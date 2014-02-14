# --- Network mesh retrieval ---

import urllib
class DataURLOpener(urllib.FancyURLopener):
    def __init__(self, url, filename):
        urllib.FancyURLopener.__init__(self)
        self.url = url
        self.filename = filename

    def retrieve(self, reporter=None, data=None):
        urllib.FancyURLopener.retrieve(self, self.url, self.filename, reporter, data)

    def http_error_default(self, url, fp, errcode, errmsg, headers):
        raise IOError(str(errcode)+" "+errmsg+", "+self.url)

def retrieve(filename, urlbase='http://simula.no/~jobh/cbcflow'):
    if not filename.endswith(".gz"):
        # Enforcing .gz extension is a quick fix to avoid trouble when
        # httpserver serves .gz file without extension, which is then
        # unreadable for dolfin.
        filename += ".gz"
    if on_master_process() and not os.path.exists(filename):
        url = urlbase+'/'+filename
        warning('%s not found, fetching from %s'%(filename,url))

        targetdir = os.path.abspath(filename[:filename.rfind('/')])
        log_level = get_log_level()
        set_log_level(PROGRESS)
        progress = [Progress(filename.split('/')[-1])]
        def reporter(numblocks, blocksize, totalsize):
            progress[0] += numblocks*blocksize / totalsize

        if not os.path.isdir(targetdir):
            os.makedirs(targetdir)
        try:
            DataURLOpener(url, filename).retrieve(reporter)
        except:
            if os.path.exists(filename):
                os.remove(filename)
            raise

        del progress[0]
        set_log_level(log_level)

    MPI.barrier()
    return filename