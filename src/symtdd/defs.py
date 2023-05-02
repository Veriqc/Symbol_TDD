from enum import Enum, auto


class DDType(Enum):
    TS = auto() # Tensornetwork with tensors
    TDD = auto()
    STDD = auto()
    TrTDD = auto()
    ExpTDD = auto()