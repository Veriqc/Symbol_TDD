from __future__ import annotations

import sys
from typing import Any

import numpy as np

if sys.version_info >= (3, 11):
    from typing import Self
else:
    from typing_extensions import Self


class Index:
    def __init__(self, *args, idx=0, hypridx=0) -> None:
        self.key = tuple(args + [idx])
        self.hypridx = hypridx

    def __eq__(self, other) -> Any:
        return self.key == other.key  # and self.idx == other.idx

    def __hash__(self) -> int:
        return hash(self.key)

    def __repr__(self) -> str:
        """ The hyper index string may exists more than two times. """
        return "_".join(str(x) for x in self.key)

    def __str__(self) -> str:
        """ The index exists no more than two times. """
        return repr(self) + "#" + str(self.hypridx)
    
    def update(self, *args, idx=0, hypridx=0) -> None:
        self.key = tuple(args + [idx])
        self.hypridx = hypridx

    def create_next(self) -> Self:
        """ Create a new Index increasing by 1 from the current. """
        type(self)(self.key[:-1], self.key[-1] + 1, self.hypridx)

    def create_next_hypr(self) -> Self:
        """ Create a new Index with hypridx increasing by 1 from the current. """
        type(self)(self.key[:-1], idx=self.key[-1], hypridx=self.hypridx + 1)


class Tensor:
    def __init__(self, data=[], indices=(), name=None) -> None:
        self.data = data
        self.indices = indices
        self.name = name

    @property
    def str_indices(self) -> list[str]:
        return [str(index) for index in self.indices]

    @property
    def str_hyprindices(self) -> list[str]:
        return [repr(index) for index in self.indices]

    def contract(self, other:Self) -> Self:
        self_indices_set = set(self.indices)
        other_indices_set = set(other.indices)
        union_indices = self_indices_set.union(other_indices_set)
        out_indices = self_indices_set.symmetric_difference(other_indices_set)

        idx_2_int = {v: i for i, v in enumerate(union_indices)}
        out_int_indices = tuple(map(idx_2_int.get, out_indices))

        def get_data_int_pair(x):
            return (x.data, tuple(map(idx_2_int.get, x.indices)))

        new_data = np.einsum(
            *get_data_int_pair(self), *get_data_int_pair(other), out_int_indices
        )
        return type(self)(new_data, tuple(out_indices))

    def tdd(self) -> None:
        pass