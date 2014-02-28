
import sys, os, glob

from commands import getstatusoutput

dry_run = False

def call(cmd):
    print "CALL:", cmd
    
    if dry_run:
        return 0, ""
    return getstatusoutput(cmd)

def call_in(path, cmd):
    print "CALL IN:", path, cmd
    if dry_run:
        return 0, ""
    p = os.path.abspath(os.curdir)
    os.chdir(path)
    try:
        s, o = getstatusoutput(cmd)
    finally:
        os.chdir(p)
    return s, o

class ReferencesRepository:
    def __init__(self, output_dir, data_dir, data_id_file, data_repo_git, data_repo_https, exclude_patterns):
        self.output_dir = output_dir
        self.data_dir = data_dir
        self.data_id_file = data_id_file
        self.data_repo_git = data_repo_git
        self.data_repo_https = data_repo_https
        self.exclude_patterns = exclude_patterns + [".git"]

    def info(self, msg):
        print msg

    def read_data_id(self):
        data_id = open(self.data_id_file, 'r').readlines()[0].strip()
        self.info("Read reference data commit id '%s' from file." % data_id)
        return data_id

    def write_data_id(self, data_id):
        self.info("Writing new reference data commit id to file.")
        open(self.data_id_file, 'w').write(str(data_id))
        s, o = call('git commit %s -m"Update reference data pointer to %s."' % (self.data_id_file, data_id))
        return s

    def clone(self):
        self.info("Executing initial cloning of reference data repository.")

        # Assuming this is only run initially
        assert not os.path.isdir(self.data_dir)

        # Try with ssh keys first
        s, o = call("git clone %s" % self.data_repo_git)

        # Fall back to https with password
        if not os.path.isdir(self.data_dir):
            s, o = call("git clone %s" % self.data_repo_https)
        return s

    def checkout(self):
        self.info("Pulling new data into existing reference data repository.")

        # Check out the master branch
        s = call_in(self.data_dir, "git checkout master")
        if s:
            self.info("Failed to checkout master, check state of reference data directory.")
            return s

        # Fetch latest data
        s, o = call_in(self.data_dir, "git fetch")
        if s:
           self.info("WARNING: Failed to fetch latest reference data from server.")
        else:
            s, o = call_in(self.data_dir, "git pull")
            if s:
                self.info("Failed to pull latest reference data from server, possibly a merge situation.")
                return 1
        return 0

    def clone_or_checkout(self):
        # Clone if it's not there, checkout if it is
        if not os.path.isdir(self.data_dir):
            s = self.clone()
        else:
            s = self.checkout()

        # Fail if it didn't work
        if s or not os.path.isdir(self.data_dir):
            self.info("Failed to get or update reference data directory '%s'." % self.data_dir)
            return 1

        return s

    def checkout_data(self, data_id=None):
        # If no commit id for the data is given, read the last one from file
        data_id = data_id or self.read_data_id()

        self.info("Trying to get data with id %s" % (data_id,))

        # If no data repository exists, clone it
        if not os.path.isdir(self.data_dir):
            s = self.clone()
            if s:
                # No repository? Nothing we can do then.
                return s

        # Try to check out the given data id without fetch from the network
        # (which takes time and may not be available)
        cmd = "git checkout -B auto %s" % (data_id,)
        s, o = call_in(self.data_dir, cmd)

        # If that failed, first update the repository (involves network access) and then try again
        if s:
            s = self.checkout()
            if not s:
                s, o = call_in(self.data_dir, cmd)

        # If it still failed, inform the user
        if s:
            self.info("Failed to checkout data with id %s." % data_id)
        return s

    def rsync(self, delete):
        self.info("Copying new reference data to %s" % self.data_dir)
        excludes = " ".join("--exclude='%s'" % e for e in self.exclude_patterns)
        flags = "--delete" if delete else ""
        cmd = "rsync -r %s %s %s/ %s" % (flags, excludes, self.output_dir, self.data_dir)
        s, o = call(cmd)
        return s

    def commit_data(self, repo_id):
        s, o = call_in(self.data_dir, "git add *")
        if not s:
            s, o = call_in(self.data_dir, 'git commit -a -m "Update reference data, current parent project head is %s."' % (repo_id,))
            self.info(o)
        return s

    def commit_data_id(self):
        s, data_id = call_in(self.data_dir, "git rev-list --max-count 1 HEAD")
        data_id = data_id.strip()
        self.info("Got data commit id '%s' from git." % data_id)
        if not s:
            s = self.write_data_id(data_id)
        return s

    def upload(self, delete):

        if not os.path.isdir(self.output_dir):
            self.info("Missing data directory '%s'." % self.output_dir)
            return 1

        # Get current id for main repo (does not include dirty files, so not quite trustworthy!)
        s, repo_id = call("git rev-list --max-count 1 HEAD")
        repo_id = repo_id.strip()
        self.info("Got repo commit id '%s' from git." % repo_id)

        # Check out latest master branch in regression repo
        s, o = call_in(self.data_dir, "git checkout master")
        if s:
            self.info("Checkout of reference repository master branch failed!")
            return s
        s, o = call_in(self.data_dir, "git pull origin master")
        if s:
            self.info("Pull of reference repository master branch failed!")
            return s

        # Copy references
        s = self.rsync(delete)
        if s:
            self.info("rsync failed!")
            return s

        # Add new files to reference data repository
        s = self.commit_data(repo_id)
        if s:
            self.info("Failed to add new files to reference data repository!")
            return s

        # Commit reference data commit id to file in main repo
        s = self.commit_data_id()
        if s:
            self.info("Failed to add data commit id to main repository!")
            return s
        
        # Push references to server
        s, o = call_in(self.data_dir, "git push origin master")
        if s:
            self.info("WARNING: Failed to push new reference data to server. You may have lost a race condition with another developer. Manual fix is necessary.")
        return s
        

