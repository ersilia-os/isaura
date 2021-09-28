class Writer(IsauraBase):

    def __init__(self, model_id, api_name):
        IsauraBase.__init__(self, model_id)
        self.model = ErsiliaModel(model_id)
        self.model.serve()

    def writer(self, api_name, input)
        pass
