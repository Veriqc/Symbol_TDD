from __future__ import annotations

import math
import sys
from typing import Any

import numpy as np

if sys.version_info >= (3, 11):
    from typing import Self
else:
    from typing_extensions import Self


class Index:
    def __init__(self, *args, hypridx=0) -> None:
        self.key = tuple(args)
        self.idx = hypridx

    def __eq__(self, other) -> Any:
        return self.key == other.key  # and self.idx == other.idx

    def __hash__(self) -> int:
        return hash(self.key)

    def __repr__(self) -> str:
        return "_".join(str(x) for x in self.key)

    def __str__(self) -> str:
        return repr(self) + "#" + str(self.idx)


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

    def contract(self, other) -> Self:
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


class TensorNetwork:
    def __init__(self, tensors=[]) -> None:
        self.tensors = tensors

    def contract(self, optimizer=None) -> Any:
        assert len(self.tensors) > 1

        if optimizer:
            raise AssertionError("Not Implemented!")
        else:
            # Direct contraction (forward)
            n = len(self.tensors)
            path = [(0, 1)] + [(0, k) for k in range(n - 2, 0, -1)]

        return self.contract_by_path(path)

    def contract_by_path(self, path) -> Any:
        self.reset_record()

        tensors = self.tensors

        for pos_pair in path:
            ts_pair = tuple(tensors[pos] for pos in pos_pair)
            new_ts = ts_pair[0].contract(ts_pair[1])

            self.record(new_ts)

            # TODO: use heapq to do in O(logn)?
            # https://stackoverflow.com/a/10163422
            tensors = [ts for i, ts in enumerate(tensors) if i not in pos_pair]
            tensors.append(new_ts)

        assert len(tensors) == 1
        return tensors[0]

    def reset_record(self) -> None:
        self.max_size = 0

    def record(self, tensor) -> None:
        ts_size = math.prod(tensor.data.shape)
        self.max_size = max(ts_size, self.max_size)
