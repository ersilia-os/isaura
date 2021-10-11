from ..core.base import IsauraBase


class Retyper(IsauraBase):
    def __init__(self, model_id):
        IsauraBase.__init__(self, model_id=model_id)

    def retype(self):
        pass
