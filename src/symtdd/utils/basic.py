import operator
from functools import reduce

from typing import Iterable, Any

def reduce_add(seq: Iterable[Any]) -> Any:
    return reduce(operator.add, seq)