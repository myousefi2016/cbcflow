
from regression_repository import ReferencesRepository

class CbcflowReferencesRepository(ReferencesRepository):
    def __init__(self):
        ReferencesRepository.__init__(self,
            output_dir = "output",
            data_dir = "cbcflow-reference-data",
            data_id_file = "cbcflow-reference-data-id",
            data_repo_git = "git@bitbucket.org:simula_cbc/cbcflow-reference-data.git",
            data_repo_https = "https://bitbucket.org/simula_cbc/cbcflow-reference-data.git",
            exclude_patterns = ["README.rst"],
            )

def download():
    "Download data from reference data repository to match data id file."
    repo = CbcflowReferencesRepository()
    repo.checkout_data()

def upload_incremental():
    """Upload data from output/ to reference data repository.

    KEEPING files from the reference data repository that are not in output/.
    """
    repo = CbcflowReferencesRepository()
    repo.upload(delete=False)

def upload_replace():
    """Upload data from output/ to reference data repository.

    DELETING files from the reference data repository that are not in output/.
    """
    repo = CbcflowReferencesRepository()
    repo.upload(delete=True)
