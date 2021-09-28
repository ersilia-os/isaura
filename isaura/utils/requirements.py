import importlib


class ErsiliaRequirement(object):

    def __init__(self):
        self.name = "ersilia"

    def import(self):
        importlib.import_module(self.name)
