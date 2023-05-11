VERSION = "0.0.1"


import sys
if sys.version_info >= (3, 11):
    from typing import Self
else:
    from typing_extensions import Self
    
import logging
from .defs import IndexType, TensorType, TDDWeightType
from .ts import Tensor, Index, HyperIndex
from .tn import TensorNetwork
from .tdd import TDD

logging.basicConfig()