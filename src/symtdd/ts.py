from __future__ import annotations

import logging
from typing import Any
import math
from collections import Counter

import numpy as np
# import opt_einsum as oe

from . import Self
from .base import IndexBase


log = logging.getLogger(__name__)
log.setLevel("DEBUG")


class Index(IndexBase):
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
    def __init__(self, data=(), indices=(), name=None) -> None:
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

    @property
    def size(self) -> int:
        return math.prod(self.data.shape)

    @property
    def index_set(self) -> Counter:
        return set(self.indices)

    @property
    def index_counter(self) -> Counter:
        return Counter(self.indices)
    
    @classmethod
    def contract_inner(cls, ts1, ts2, out_indices, intersect_indices, union_indices):
        # out_int_indices = out_indices
        # We still use idx_2_int to avoid number limit of indice from either numpy or opt_einsum
        
        idx_2_int = {v: i for i, v in enumerate(union_indices)}
        out_int_indices = tuple(map(idx_2_int.get, out_indices))

        def get_data_ints_pair(x):
            y = tuple(map(idx_2_int.get, x.indices))
            log.debug("\ninput: %s\noutput: %s", x, y)
            return (x.data, y)
        
        # self_di_pair = (self.data, self.indices)
        # other_di_pair = (other.data, other.indices)
        self_di_pair = get_data_ints_pair(ts1)
        other_di_pair = get_data_ints_pair(ts2)

        # new_data = oe.contract(
        #     *self_di_pair, *other_di_pair, out_int_indices
        # )
        new_data = np.einsum(
            *self_di_pair, *other_di_pair, out_int_indices
        )

        # log.debug("data: %s", new_data)
        log.debug("data shape: %s", new_data.shape)

        return cls(new_data, tuple(out_indices))

    def contract(self, other: Self, index_counter: Counter | dict[HyperIndex, int]=None) -> Self:
        self_indices_set = self.index_set
        other_indices_set = other.index_set

        intersect_indices = self_indices_set.intersection(other_indices_set)
        union_indices = self_indices_set.union(other_indices_set)

        out_indices = self_indices_set.symmetric_difference(other_indices_set) # = union - intersect

        if index_counter is not None:
            # Handle hyperindex contraction
            # log.debug("global index counter: %s", index_counter)
            self_counter = self.index_counter
            other_counter = other.index_counter
            for index in intersect_indices:
                count = self_counter[index] + other_counter[index]
                if count != index_counter[index]:
                    out_indices.add(index)
                    count -= 1
                index_counter[index] -= count

        # Warn: out_indices is a set and is in arbitary order!
        log.debug("====Contract")
        log.debug("%s,%s->%s", self.indices, other.indices, out_indices)

        tensor = self.contract_inner(self, other, out_indices, intersect_indices, union_indices)

        log.debug("Contract End====")

        return tensor
    
    def sort(self) -> None:
        sort_idxs = np.argsort(self.indices)
        self.data = np.moveaxis(self.data, sort_idxs, np.arange(self.data.ndim))
        self.indices = tuple(self.indices[idx] for idx in sort_idxs)

    def tdd(self) -> None:
        pass