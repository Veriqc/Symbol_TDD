from enum import Enum, auto


class IndexType(Enum):
    STANDARD = auto()
    HYPEREDGE = auto()


class TensorType(Enum):
    TENSOR = auto()
    TDD = auto()


class TDDWeightType(Enum):
    NA = auto()
    COMPLEX = auto()
    SymDD = auto()
    TrDD = auto()
    ExpDD = auto()