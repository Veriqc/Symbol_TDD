from __future__ import annotations

from . import Self, Tensor, TDDWeightType
from .base import DDBase, NodeBase, EdgeBase

import numpy as np


class Node(NodeBase):
    pass


class Edge(EdgeBase):
    pass


class TDD(Tensor, DDBase):
    unique_table = {}
    compute_table = {}
    terminal = Node.terminal()

    def __init__(self, tensor: Tensor, edge: Edge) -> None:
        # self.data = tensor.data
        self.root_edge = edge
        # TODO: unify key_2_idx (dict) or indices (list)
        self.indices = tensor.indices
        self.name = tensor.name

    def __str__(self) -> str:
        return f"TDD {self.name}: indices {self.str_hyprindices}"

    @property
    def key_2_idx(self):
        return self.indices
    
    @classmethod
    def contract_inner(cls, tdd1, tdd2, out_indices, union_indices, intersect_indices):
        raise ValueError("Not implemented")

    @classmethod
    def from_tensor(cls, tensor: Tensor, wtype: TDDWeightType) -> Self:
        if wtype is not TDDWeightType.COMPLEX:
            raise ValueError("Not implemented")
        
        # Hack! Direct access tensor member
        def np_2_tdd_recur(data: np.ndarray | complex, d: int) -> Edge:
            if not isinstance(data, np.ndarray):
                return Edge(data, cls.terminal)
            else:
                size = data.shape[0]
                succ = [np_2_tdd_recur(data[i], d + 1) for i in range(size)]
                return Edge(1.0+0j, Node(d, succ))
                
        edge = np_2_tdd_recur(tensor.data, 0)
        
        return cls(tensor, edge)