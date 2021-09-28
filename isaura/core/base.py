import os
from ..default import REPOSITORY_PATH
from ..default import HDF5_EXTENSION
from .. import logger


class IsauraBase(object):

    def __init__(self, model_id):
        self.logger = logger
        self.model_id = model_id
        self.repository_path = REPOSITORY_PATH

    @property
    def data_path(self):
        return os.path.join(repository_path, self.model_id+HDF5_EXTENSION)
