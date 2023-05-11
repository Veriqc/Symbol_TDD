from __future__ import annotations

from graphviz import Digraph
from IPython.display import Image

from typing import Iterable

from . import Self


class IndexBase():
    pass


class IndirectDigraph(Digraph):
    def __init__(self, key_2_idx, **kargs) -> None:
        super().__init__(**kargs)
        self.key_2_idx = key_2_idx


class NodeBase():
    CONFIG = dict(fontname="helvetica", shape="circle", color="red")

    def __init__(self, key: int, successors: Iterable[EdgeBase]=()) -> None:
        self.key = key
        self.successors = successors

    @classmethod
    def terminal(cls) -> Self:
        return cls(-1)

    @property
    def is_terminal(self) -> bool:
        return self.key == -1

    @property
    def name(self) -> str:
        return str(id(self))
    
    def get_label(self, key_2_idx: dict[int, IndexBase]) -> str:
        return str(key_2_idx[self.key])

    def draw(self, dot: IndirectDigraph) -> None:
        if self.is_terminal:
            dot.node(str(-1), str(1), **self.CONFIG)
        else:
            dot.node(self.name, self.get_label(dot.key_2_idx), **self.CONFIG)


class EdgeBase():
    def __init__(self, weight: complex, node: NodeBase) -> None:
        self.w = weight
        self.v = node

    def draw(self, dot: IndirectDigraph, src_name: str) -> None:
        dot.edge(src_name, self.v.name, color="blue", label=str(self.w))


class DDBase():
    def __init__(self, root_edge: EdgeBase, key_2_idx: dict[int, IndexBase]) -> None:
        self.root_edge = root_edge
        self.key_2_idx = key_2_idx

    def show(self, real_label=True, filename="output"):
        dot = IndirectDigraph(name="reduced_tree")
        dot.format = "png"
        dot.node("o", "", shape="none")
        layout(dot, "o", self.root_edge)
        return Image(dot.render(filename))
    

def layout(dot: IndirectDigraph, src_name: str, edge: EdgeBase) -> None:
    edge.draw(dot, src_name)
    v = edge.v
    v.draw(dot)
    for succ in v.successors:
        layout(dot, v.name, succ)