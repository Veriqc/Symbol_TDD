from __future__ import annotations

from . import Self, Tensor, TDDWeightType


class TDD(Tensor):
    @classmethod
    def from_tensor(tensor: Tensor, wtype: TDDWeightType) -> Self:
        return tensor