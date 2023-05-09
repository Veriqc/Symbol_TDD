from __future__ import annotations

import sys
import logging
from typing import Any
from collections import Counter

import numpy as np
# import opt_einsum as oe

if sys.version_info >= (3, 11):
    from typing import Self
else:
    from typing_extensions import Self


log = logging.getLogger(__name__)
log.setLevel("DEBUG")


class Index:
    """
    Though here we implement hyper indices storage,
    the indices we compare (__eq__, __lt__) and represent(__repr__, __str__)
    should still be with standard. (Each index exists no more than two times.)
    """
    def __init__(self, *args, idx=0, hypridx=0) -> None:
        self.key = args + (idx,)
        self.hypridx = hypridx

    def __eq__(self, other) -> Any:
        return self.key == other.key and self.hypridx == other.hypridx
    
    def __lt__(self, other) -> bool:
        for a, b in zip(self.key, other.key):
            if a < b:
                return True
            elif a > b:
                return False
        return self.hypridx < other.hypridx

    def __hash__(self) -> int:
        return hash(self.key)

    def __repr__(self) -> str:
        return self._str + "#" + str(self.hypridx)

    def __str__(self) -> str:
        return self._str + "_" + str(self.hypridx)
    
    @property
    def _str(self) -> str:
        return "_".join(str(x) for x in self.key)
    
    def update(self, *args, idx=0, hypridx=0) -> None:
        self.key = args + (idx,)
        self.hypridx = hypridx

    def create_next(self, with_hypridx=False) -> Self:
        """ Create a new Index increasing by 1 from the current. """
        pre_args = self.key[:-1]
        idx, hypridx = self.key[-1], self.hypridx
        if with_hypridx:
            hypridx += 1
        else:
            idx += 1
        new_index = type(self)(*pre_args, idx=idx, hypridx=hypridx)
        return new_index
    

class HyperIndex(Index):
    __hash__ = Index.__hash__

    def __init__(self, index:Index) -> None:
        self.key = index.key
        self.hypridx = index.hypridx

    def __eq__(self, other) -> Any:
        return self.key == other.key

    def __repr__(self) -> str:
        return self._str

    def __str__(self) -> str:
        return self._str


class Tensor:
    def __init__(self, data=[], indices=(), name=None) -> None:
        self.data = data
        self.indices = indices
        self.name = name

    def __str__(self) -> str:
        return f"TS {self.name}: shape {self.data.shape}, indices {self.str_hyprindices}"

    @property
    def str_indices(self) -> list[str]:
        return [str(index) for index in self.indices]

    @property
    def str_hyprindices(self) -> list[str]:
        return [repr(index) for index in self.indices]

    def contract(self, other:Self, order_counts: list[int]) -> Self:
        log.debug("order counter: %s", order_counts)
        # Since we use indirect index, indice here are integers (global orders).

        self_counter = Counter(self.indices)
        other_counter = Counter(other.indices)

        self_indices_set = set(self.indices)
        other_indices_set = set(other.indices)
        intersect_indices = self_indices_set.intersection(other_indices_set)
        out_indices = self_indices_set.symmetric_difference(other_indices_set)
        for index in intersect_indices:
            count = self_counter[index] + other_counter[index]
            if count != order_counts[index]:
                out_indices.add(index)
                count -= 1
            order_counts[index] -= count

        # out_int_indices = out_indices
        # We still use idx_2_int to avoid number limit of indice from either numpy or opt_einsum

        union_indices = self_indices_set.union(other_indices_set)
        idx_2_int = {v: i for i, v in enumerate(union_indices)}
        out_int_indices = tuple(map(idx_2_int.get, out_indices))

        # log.debug("====Contract")
        log.debug("====Contract\nmap: %s", idx_2_int)

        def get_data_ints_pair(x):
            y = tuple(map(idx_2_int.get, x.indices))
            log.debug("\ninput: %s\noutput: %s", x, y)
            return (x.data, y)
        
        # self_di_pair = (self.data, self.indices)
        # other_di_pair = (other.data, other.indices)
        self_di_pair = get_data_ints_pair(self)
        other_di_pair = get_data_ints_pair(other)

        # new_data = oe.contract(
        #     *self_di_pair, *other_di_pair, out_int_indices
        # )
        new_data = np.einsum(
            *self_di_pair, *other_di_pair, out_int_indices
        )

        log.debug("%s,%s->%s", self_di_pair[1], other_di_pair[1], out_int_indices)
        # log.debug("data: %s", new_data)
        log.debug("Contract End====")

        return type(self)(new_data, tuple(out_indices))
    
    def sort(self) -> None:
        sort_idxs = np.argsort(self.indices)
        self.data = np.moveaxis(self.data, sort_idxs, np.arange(self.data.ndim))
        self.indices = tuple(self.indices[idx] for idx in sort_idxs)

    def tdd(self) -> None:
        pass