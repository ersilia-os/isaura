from .base import IsauraBase


class IsauraApi(IsauraBase):
    def __init__(self, model_id):
        self.model_id = model_id

    def get(self, inputs):
        pass
