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

    def __init__(self, edge: Edge, indices=(), name='') -> None:
        # self.data = tensor.data
        self.root_edge = edge
        # TODO: unify key_2_idx (dict) or indices (list)
        self.indices = indices
        self.name = name

    def __str__(self) -> str:
        return f"TDD {self.name}: indices {self.str_hyprindices}"

    @property
    def key_2_idx(self):
        return self.indices
    
    @property
    def is_trivial(self) -> bool:
        return self.root_edge.is_terminated
    
    @classmethod
    def contract_inner(cls, f: TDD, g: TDD, out_indices, intersect_indices):
        raise ValueError("Not implemented yet")
        if len(intersect_indices) == 0:
            return cls.mul(f, g, out_indices)
        
        x = intersect_indices[0]
        f0 = f.slice(x, 0)
        f1 = f.slice(x, 1)
        g0 = g.slice(x, 0)
        g1 = g.slice(x, 1)
        l = cls.contract_inner(f0, g0, out_indices, union_indices, intersect_indices[1:])
        r = cls.contract_inner(f1, g1, out_indices, union_indices, intersect_indices[1:])
        return cls.add(l, r, out_indices)


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
        
        return cls(edge, tensor.indices, tensor.name)