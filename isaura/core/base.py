import os
from ..default import REPOSITORY_PATH
from ..default import HDF5_EXTENSION
from .. import logger


class IsauraBase(object):

    def __init__(self, model_id, verbose=True):
        self.logger = logger
        if verbose:
            self.logger.set_verbosity(1)
        else:
            self.logger.set_verbosity(0)
        self.model_id = model_id
        self.repository_path = REPOSITORY_PATH

    @property
    def file_name(self):
        return self.model_id + "." + HDF5_EXTENSION

    @property
    def data_path(self):
        return os.path.join(self.repository_path, self.file_name)

    def exists(self):
        return os.path.exists(self.data_path)
