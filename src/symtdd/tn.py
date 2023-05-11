from __future__ import annotations

import sys
from typing import Any

if sys.version_info >= (3, 11):
    from typing import Self
else:
    from typing_extensions import Self

from .ts import Tensor
from .utils.basic import reduce_add


class TensorNetwork:
    def __init__(self, tensors=(), usehyper=False) -> None:
        """
        usehyper: determine if tn should contract by hyperindice (each index exists more than two times)
        """
        self.tensors = tensors
        self.usehyper = usehyper

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
        index_counter = reduce_add([tensor.index_counter for tensor in tensors]) if self.usehyper else None

        for pos_pair in path:
            ts_pair = tuple(tensors[pos] for pos in pos_pair)
            new_ts = ts_pair[0].contract(ts_pair[1], index_counter)

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
        self.max_size = max(tensor.size, self.max_size)
