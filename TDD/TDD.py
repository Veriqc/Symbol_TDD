"""Dummy module docstring"""
import copy

import numpy as np
from graphviz import Digraph
from IPython.display import Image
from sympy import *

# Define global variables
computed_table = dict()
unique_table = dict()
global_index_order = dict()
global_node_idx = 0
add_find_time = 0
add_hit_time = 0
cont_find_time = 0
cont_hit_time = 0
epi = 1e-10


class Index:
    """The index, here idx is used when there is a hyperedge"""

    def __init__(self, key, idx=0):
        self.key = key
        self.idx = idx

    def __eq__(self, other):
        if self.key == other.key and self.idx == other.idx:
            return True
        else:
            return False

    def __lt__(self, other):
        if global_index_order[self.key] < global_index_order[other.key]:
            return True
        elif self.key == other.key and self.idx < other.idx:
            return True

        return False

    def __str__(self):
        return str((self.key, self.idx))


class Node:
    """To define the node of TDD"""

    def __init__(self, key, num=2):
        self.idx = 0
        self.key = key
        self.succ_num = num
        self.out_weight = []
        for k in range(num):
            self.out_weight.append(S_one)
        self.successor = [None] * num
        self.meas_prob = []

    def __repr__(self) -> str:
        return str(self.key) + str(self.out_weight) + str(self.successor)

    def mystr(self, depth=0) -> str:
        return (
            "├───→" * depth
            + "key:"
            + str(self.key)
            + " out_weight:"
            + str(self.out_weight)
            + "\n "
            + "".join(
                [
                    self.successor[i].mystr(depth + 1)
                    if self.successor[i] is not None
                    else ""
                    for i in range(self.succ_num)
                ]
            )
        )


class TDD:
    def __init__(self, node):
        """TDD"""
        self.weight = S_one

        self.index_set = []

        self.key_2_index = dict()
        self.index_2_key = dict()
        self.key_width = dict()  # Only used when change TDD to numpy

        if isinstance(node, Node):
            self.node = node
        else:
            self.node = Node(node)

    def mystr(self) -> str:
        return str(self.weight) + "\n" + self.node.mystr()

    def node_number(self):
        node_set = set()
        node_set = get_node_set(self.node, node_set)
        return len(node_set)

    def self_copy(self):
        temp = TDD(self.node)
        temp.weight = self.weight
        temp.index_set = copy.copy(self.index_set)
        temp.key_2_index = copy.copy(self.key_2_index)
        temp.index_2_key = copy.copy(self.index_2_key)
        return temp

    def show(self, real_label=True, filename="output"):
        edge = []
        dot = Digraph(name="reduced_tree")
        dot = layout(self.node, self.key_2_index, dot, edge, real_label)
        dot.node("o", "", shape="none")
        label1 = str(self.weight)
        dot.edge("o", str(self.node.idx), color="blue", label=label1)
        dot.format = "png"
        return Image(dot.render(filename))

    def __eq__(self, other):
        return self.node == other.node and self.weight == other.weight


def layout(node, key_2_idx, dot=Digraph(), succ=[], real_label=True):
    col = ["red", "blue", "black", "green"]
    if real_label and node.key in key_2_idx:
        if node.key == -1:
            dot.node(
                str(node.idx), str(1), fontname="helvetica", shape="circle", color="red"
            )
        else:
            dot.node(
                str(node.idx),
                key_2_idx[node.key],
                fontname="helvetica",
                shape="circle",
                color="red",
            )
    else:
        dot.node(
            str(node.idx),
            str(node.key),
            fontname="helvetica",
            shape="circle",
            color="red",
        )
    for k in range(node.succ_num):
        if node.successor[k]:
            #             label1=str(node.out_weight[k])
            label1 = str(node.out_weight[k])
            if not node.successor[k] in succ:
                dot = layout(node.successor[k], key_2_idx, dot, succ, real_label)
                #                 if get_int_key(node.out_weight[k].weight)==(0,0):
                #                     continue
                dot.edge(
                    str(node.idx),
                    str(node.successor[k].idx),
                    color=col[k % 4],
                    label=label1,
                )
                succ.append(node.successor[k])
            else:
                #                 if get_int_key(node.out_weight[k].weight)==(0,0):
                #                     continue
                dot.edge(
                    str(node.idx),
                    str(node.successor[k].idx),
                    color=col[k % 4],
                    label=label1,
                )
    return dot


def Ini_TDD(index_order=[], var=[], n=50, type=None, unique_table_reset=True):
    """To initialize the unique_table,computed_table and set up a global index order"""
    global computed_table
    global unique_table
    global global_node_idx
    global add_find_time, add_hit_time, cont_find_time, cont_hit_time
    global S_one, S_zero
    global tdd_type

    if type == "SymTDD" or type == None:
        from TDD.SymTDD.BDD import Ini_BDD, get_bdd

        tdd_type = "SymTDD"
    if type == "TrDD":
        from TDD.TrDD.BDD import Ini_BDD, get_bdd

        tdd_type = "TrDD"
    if type == "Exp":
        from TDD.Exp.EXP import Ini_BDD, get_bdd

        tdd_type = "Exp"

    if unique_table_reset == True:
        global_node_idx = 0
        unique_table = dict()
        computed_table = dict()

        if tdd_type:
            if not var:
                var = ["x" + str(k) for k in range(n - 1, -1, -1)]
            Ini_BDD(var)
    add_find_time = 0
    add_hit_time = 0
    cont_find_time = 0
    cont_hit_time = 0
    set_index_order(index_order)

    S_one = get_bdd(1)
    S_zero = get_bdd(0)

    return get_identity_tdd()


def Clear_TDD():
    """To initialize the unique_table,computed_table and set up a global index order"""
    global computed_table
    global unique_table
    global global_node_idx
    global add_find_time, add_hit_time, cont_find_time, cont_hit_time
    global_node_idx = 0
    unique_table.clear()
    computed_table.clear()
    add_find_time = 0
    add_hit_time = 0
    cont_find_time = 0
    cont_hit_time = 0
    global_node_idx = 0


#     return 1


def get_identity_tdd():
    node = Find_Or_Add_Unique_table(-1)
    tdd = TDD(node)
    tdd.index_2_key = {-1: -1}
    tdd.key_2_index = {-1: -1}
    return tdd


def get_zero_tdd():
    node = Find_Or_Add_Unique_table(-1)
    tdd = TDD(node)
    tdd.weight = S_zero
    tdd.index_2_key = {-1: -1}
    tdd.key_2_index = {-1: -1}
    return tdd


def get_unique_table():
    return unique_table


def get_unique_table_num():
    return len(unique_table)


def set_index_order(var_order):
    global global_index_order
    global_index_order = dict()
    if isinstance(var_order, list):
        for k in range(len(var_order)):
            global_index_order[var_order[k]] = k
    if isinstance(var_order, dict):
        global_index_order = copy.copy(var_order)
    global_index_order[-1] = float("inf")


def get_index_order():
    global global_index_order
    return copy.copy(global_index_order)


def get_int_key(weight):
    """To transform a complex number to a tuple with int values"""
    global epi

    return (int(round(weight.real / epi)), int(round(weight.imag / epi)))


def get_node_set(node, node_set=set()):
    """Only been used when counting the node number of a TDD"""
    if not node in node_set:
        node_set.add(node)
        for k in range(node.succ_num):
            if node.successor[k]:
                node_set = get_node_set(node.successor[k], node_set)
    return node_set


def get_weight(sl):
    return (get_int_key(sl.weight), sl.node)


def Find_Or_Add_Unique_table(x, weigs=[], succ_nodes=[]):
    """To return a node if it already exist, creates a new node otherwise"""
    global global_node_idx, unique_table

    if x == -1:
        if unique_table.__contains__(x):
            return unique_table[x]
        else:
            res = Node(x)
            res.idx = 0
            unique_table[x] = res
        return res
    temp_key = [x]
    for k in range(len(weigs)):
        temp_key.append(get_weight(weigs[k]))
        temp_key.append(succ_nodes[k])

    temp_key = tuple(temp_key)
    if temp_key in unique_table:
        return unique_table[temp_key]
    else:
        res = Node(x, len(succ_nodes))
        global_node_idx += 1
        res.idx = global_node_idx
        res.out_weight = weigs
        res.successor = succ_nodes
        unique_table[temp_key] = res
    return res


def if_line_combine(tdd1, tdd2):
    if tdd1.node == tdd2.node:
        return False
    if tdd1.node.key != tdd2.node.key:
        return False
    if (
        tdd1.node.successor[0]
        == tdd1.node.successor[1]
        == tdd2.node.successor[0]
        == tdd2.node.successor[1]
    ):
        return True
    if (
        tdd1.node.successor[0] != tdd1.node.successor[1]
        or tdd2.node.successor[0] != tdd2.node.successor[1]
    ):
        return False
    if tdd1.node.out_weight != tdd2.node.out_weight:
        return False
    return if_line_combine(
        Slicing(tdd1, tdd1.node.key, 0), Slicing(tdd2, tdd2.node.key, 0)
    )


def if_line_combine2(tdd1, tdd2):
    if tdd1.node.key != tdd2.node.key:
        return False
    if tdd1.node == tdd2.node:
        return False

    if (
        tdd1.node.successor[0] == tdd2.node.successor[0]
        and tdd1.node.successor[1] == tdd2.node.successor[1]
    ):
        # print('TDD 280', tdd1.node.successor[0])
        # print('TDD 281', tdd1.node.successor[1])
        return True
    if (
        tdd1.node.successor[0] != tdd1.node.successor[1]
        and tdd1.node.successor[0].key != -1
        and tdd1.node.successor[1].key != -1
    ):
        return False
    if (
        tdd2.node.successor[0] != tdd2.node.successor[1]
        and tdd2.node.successor[0].key != -1
        and tdd2.node.successor[1].key != -1
    ):
        return False
    if tdd1.node.successor[0].key != -1:
        temp1 = Slicing(tdd1, tdd1.node.key, 0)
    else:
        temp1 = Slicing(tdd1, tdd1.node.key, 1)
    if tdd2.node.successor[0].key != -1:
        temp2 = Slicing(tdd2, tdd2.node.key, 0)
    else:
        temp2 = Slicing(tdd2, tdd2.node.key, 1)
    return if_line_combine2(temp1, temp2)


def if_line_combine3(tdd1, tdd2):
    if tdd1.node.key != tdd2.node.key:
        return False
    if tdd1.node == tdd2.node:
        return False


def normalise_line_combine(x, tdd1, tdd2):
    if (
        tdd1.node.successor[0]
        == tdd1.node.successor[1]
        == tdd2.node.successor[0]
        == tdd2.node.successor[1]
    ):
        f0 = mul_weight(tdd1.weight, tdd1.node.out_weight[0])
        f1 = mul_weight(tdd1.weight, tdd1.node.out_weight[1])
        if mul_weight(f0, f1) != S_zero:
            return normalize(x, [tdd1, tdd2], combine=False)
        g0 = mul_weight(tdd2.weight, tdd2.node.out_weight[0])
        if mul_weight(f0, g0) != S_zero:
            return normalize(x, [tdd1, tdd2], combine=False)
        if mul_weight(f1, g0) != S_zero:
            return normalize(x, [tdd1, tdd2], combine=False)
        g1 = mul_weight(tdd2.weight, tdd2.node.out_weight[1])
        if mul_weight(f0, g1) != S_zero:
            return normalize(x, [tdd1, tdd2], combine=False)
        if mul_weight(f1, g1) != S_zero:
            return normalize(x, [tdd1, tdd2], combine=False)
        if mul_weight(g0, g1) != S_zero:
            return normalize(x, [tdd1, tdd2], combine=False)
        node1 = tdd1.node.successor[0]
        node2 = tdd1.node.successor[0]
        w = S_one
        if f0 + g0 == S_zero:
            node1 = Find_Or_Add_Unique_table(-1)
        if f1 + g1 == S_zero:
            node2 = Find_Or_Add_Unique_table(-1)
        if f0 + g0 == f1 + g1:
            w = f0 + g0
            node = node1
        else:
            node = Find_Or_Add_Unique_table(
                tdd1.node.key, [f0 + g0, f1 + g1], [node1, node2]
            )
        w1 = (f0 + f1) * w
        w2 = (g0 + g1) * w
        node1 = node
        node2 = node
        w = S_one
        if w1 == S_zero:
            node1 = Find_Or_Add_Unique_table(-1)
        if w2 == S_zero:
            node2 = Find_Or_Add_Unique_table(-1)
        if w1 == w2:
            w = w1
            node = node1
        else:
            node = Find_Or_Add_Unique_table(x, [w1, w2], [node1, node2])
        tdd = TDD(node)
        tdd.weight = w
        return tdd
    temp1 = Slicing(tdd1, tdd1.node.key, 0)
    temp1.weight = tdd1.weight
    temp2 = Slicing(tdd2, tdd1.node.key, 0)
    temp2.weight = tdd2.weight
    temp_tdd = normalise_line_combine(x, temp1, temp2)
    succ0 = Slicing(temp_tdd, x, 0).node
    succ1 = Slicing(temp_tdd, x, 1).node
    temp_weight = tdd1.node.out_weight
    if temp_weight[0] == S_zero:
        succ0 = Find_Or_Add_Unique_table(-1)
    if temp_weight[1] == S_zero:
        succ1 = Find_Or_Add_Unique_table(-1)
    node1 = Find_Or_Add_Unique_table(tdd1.node.key, temp_weight, [succ0, succ0])
    node2 = Find_Or_Add_Unique_table(tdd1.node.key, temp_weight, [succ1, succ1])
    if temp_tdd.weight != S_one:
        w1 = temp_tdd.weight
        w2 = temp_tdd.weight
    else:
        w1 = temp_tdd.node.out_weight[0]
        w2 = temp_tdd.node.out_weight[1]
    if w1 == S_zero:
        node1 = Find_Or_Add_Unique_table(-1)
    if w2 == S_zero:
        node2 = Find_Or_Add_Unique_table(-1)
    w = S_one
    if node1 == node2 and w1 == w2:
        w = w1
        node = node1
    else:
        node = Find_Or_Add_Unique_table(x, [w1, w2], [node1, node2])
    tdd = TDD(node)
    tdd.weight = w
    return tdd


def normalise_line_combine2(x, tdd1, tdd2):
    #     print(x)
    #     if (tdd1.node.successor[0]==tdd2.node.successor[0] or tdd1.node.out_weight[0]==S_zero or tdd2.node.out_weight[0]==S_zero) and (tdd1.node.successor[1]==tdd2.node.successor[1] or tdd1.node.out_weight[1]==S_zero or tdd2.node.out_weight[1]==S_zero):
    #     if tdd1.node.successor==tdd2.node.successor:
    con = False
    if tdd1.node.successor == tdd2.node.successor:
        con = True
    nodes = []
    if tdd1.node.out_weight[0] != S_zero:
        nodes.append(tdd1.node.successor[0])
    if tdd1.node.out_weight[1] != S_zero:
        nodes.append(tdd1.node.successor[1])
    if tdd2.node.out_weight[0] != S_zero:
        nodes.append(tdd2.node.successor[0])
    if tdd2.node.out_weight[1] != S_zero:
        nodes.append(tdd2.node.successor[1])
    flag = True
    for k in range(len(nodes) - 1):
        if not nodes[k + 1] == nodes[0]:
            flag = False
    if con or flag:
        f0 = mul_weight(tdd1.weight, tdd1.node.out_weight[0])
        g0 = mul_weight(tdd2.weight, tdd2.node.out_weight[0])
        f1 = mul_weight(tdd1.weight, tdd1.node.out_weight[1])
        g1 = mul_weight(tdd2.weight, tdd2.node.out_weight[1])
        if mul_weight(f0 + f1, g0 + g1) != S_zero:
            return normalize(x, [tdd1, tdd2], combine=False)
        if tdd1.node.out_weight[0] != S_zero:
            node1 = tdd1.node.successor[0]
        else:
            node1 = tdd2.node.successor[0]
        if tdd1.node.out_weight[1] != S_zero:
            node2 = tdd1.node.successor[1]
        else:
            node2 = tdd2.node.successor[1]
        w = S_one
        if f0 + g0 == S_zero:
            node1 = Find_Or_Add_Unique_table(-1)
        if f1 + g1 == S_zero:
            node2 = Find_Or_Add_Unique_table(-1)
        if f0 + g0 == f1 + g1 and node1 == node2:
            w = f0 + g0
            node = node1
        else:
            node = Find_Or_Add_Unique_table(
                tdd1.node.key, [f0 + g0, f1 + g1], [node1, node2]
            )
        w1 = (f0 + f1) * w
        w2 = (g0 + g1) * w
        node1 = node
        node2 = node
        w = S_one
        if w1 == S_zero:
            node1 = Find_Or_Add_Unique_table(-1)
        if w2 == S_zero:
            node2 = Find_Or_Add_Unique_table(-1)
        if w1 == w2:
            w = w1
            node = node1
        else:
            node = Find_Or_Add_Unique_table(x, [w1, w2], [node1, node2])
        tdd = TDD(node)
        tdd.weight = w
        return tdd
    if tdd1.node.successor[0].key != -1:
        temp1 = Slicing(tdd1, tdd1.node.key, 0)
    else:
        temp1 = Slicing(tdd1, tdd1.node.key, 1)
    if tdd2.node.successor[0].key != -1:
        temp2 = Slicing(tdd2, tdd2.node.key, 0)
    else:
        temp2 = Slicing(tdd2, tdd2.node.key, 1)

    temp1.weight = tdd1.weight
    temp2.weight = tdd2.weight

    temp_tdd = normalise_line_combine2(x, temp1, temp2)

    succ0 = Slicing(temp_tdd, x, 0).node
    succ1 = Slicing(temp_tdd, x, 0).node
    temp_weight1 = tdd1.node.out_weight
    if temp_weight1[0] == S_zero:
        succ0 = Find_Or_Add_Unique_table(-1)
    if temp_weight1[1] == S_zero:
        succ1 = Find_Or_Add_Unique_table(-1)
    node1 = Find_Or_Add_Unique_table(tdd1.node.key, temp_weight1, [succ0, succ1])

    succ0 = Slicing(temp_tdd, x, 1).node
    succ1 = Slicing(temp_tdd, x, 1).node
    temp_weight2 = tdd2.node.out_weight
    if temp_weight2[0] == S_zero:
        succ0 = Find_Or_Add_Unique_table(-1)
    if temp_weight2[1] == S_zero:
        succ1 = Find_Or_Add_Unique_table(-1)
    node2 = Find_Or_Add_Unique_table(tdd2.node.key, temp_weight2, [succ0, succ1])

    #     if (node1.successor[0]==node2.successor[0] or node1.out_weight[0]==S_zero or node2.out_weight[0]==S_zero) and (node1.successor[1]==node2.successor[1] or node1.out_weight[1]==S_zero or node2.out_weight[1]==S_zero):
    con = False
    if node1.successor == node2.successor:
        con = True
    nodes = []
    if node1.out_weight[0] != S_zero:
        nodes.append(node1.successor[0])
    if node1.out_weight[1] != S_zero:
        nodes.append(node1.successor[1])
    if node2.out_weight[0] != S_zero:
        nodes.append(node2.successor[0])
    if node2.out_weight[1] != S_zero:
        nodes.append(node2.successor[1])
    flag = True
    for k in range(len(nodes) - 1):
        if not nodes[k + 1] == nodes[0]:
            flag = False
    if con or flag:
        temp1 = TDD(node1)
        temp1.weight = temp_tdd.weight * temp_tdd.node.out_weight[0]
        temp2 = TDD(node2)
        temp2.weight = temp_tdd.weight * temp_tdd.node.out_weight[1]
        return normalise_line_combine2(x, temp1, temp2)
        #         f0=mul_weight(temp1.weight,temp1.node.out_weight[0])
    #         g0=mul_weight(temp2.weight,temp2.node.out_weight[0])
    #         f1=mul_weight(temp1.weight,temp1.node.out_weight[1])
    #         g1=mul_weight(temp2.weight,temp2.node.out_weight[1])
    #         if mul_weight(f0+f1,g0+g1)!=S_zero:
    #             return normalize(x,[temp1,temp2],combine=False)
    #         if temp1.node.out_weight[0]!=S_zero:
    #             node1=temp1.node.successor[0]
    #         else:
    #             node1=temp2.node.successor[0]
    #         if temp1.node.out_weight[1]!=S_zero:
    #             node2=temp1.node.successor[1]
    #         else:
    #             node2=temp2.node.successor[1]
    #         w=S_one
    #         if f0+g0==S_zero:
    #             node1=Find_Or_Add_Unique_table(-1)
    #         if f1+g1==S_zero:
    #             node2=Find_Or_Add_Unique_table(-1)
    #         if f0+g0==f1+g1 and node1==node2:
    #             w=f0+g0
    #             node=node1
    #         else:
    #             node=Find_Or_Add_Unique_table(temp1.node.key,[f0+g0,f1+g1],[node1,node2])
    #         w1=(f0+f1)*w
    #         w2=(g0+g1)*w
    #         node1=node
    #         node2=node
    #         w=S_one
    #         if w1==S_zero:
    #             node1=Find_Or_Add_Unique_table(-1)
    #         if w2==S_zero:
    #             node2=Find_Or_Add_Unique_table(-1)
    #         if w1==w2:
    #             w=w1
    #             node=node1
    #         else:
    #             node=Find_Or_Add_Unique_table(x,[w1,w2],[node1,node2])
    #         tdd=TDD(node)
    #         tdd.weight=w
    #         return tdd

    if temp_tdd.weight != S_one:
        w1 = temp_tdd.weight
        w2 = temp_tdd.weight
    else:
        w1 = temp_tdd.node.out_weight[0]
        w2 = temp_tdd.node.out_weight[1]
    if w1 == S_zero:
        node1 = Find_Or_Add_Unique_table(-1)
    if w2 == S_zero:
        node2 = Find_Or_Add_Unique_table(-1)
    w = S_one
    if node1 == node2 and w1 == w2:
        w = w1
        node = node1
    else:
        node = Find_Or_Add_Unique_table(x, [w1, w2], [node1, node2])
    tdd = TDD(node)
    tdd.weight = w
    return tdd


def normalize(x, the_successors, combine=False):
    """The normalize and reduce procedure"""
    global epi
    all_equal = True
    for k in range(1, len(the_successors)):
        if the_successors[k] != the_successors[0]:
            all_equal = False
            break
    if all_equal:
        return the_successors[0]
    if combine:
        if if_line_combine2(the_successors[0], the_successors[1]):
            # print('TDD 557')
            return normalise_line_combine2(x, the_successors[0], the_successors[1])
    [a, b, c] = the_successors[1].weight.self_normalize(the_successors[0].weight)
    # [a,b,c]=[S_one,the_successors[0].weight,the_successors[1].weight]

    succ_nodes = [succ.node for succ in the_successors]
    node = Find_Or_Add_Unique_table(x, [b, c], succ_nodes)
    res = TDD(node)
    res.weight = a
    return res


def get_count():
    global add_find_time, add_hit_time, cont_find_time, cont_hit_time
    print("add:", add_hit_time, "/", add_find_time, "/", add_hit_time / add_find_time)
    print(
        "cont:", cont_hit_time, "/", cont_find_time, "/", cont_hit_time / cont_find_time
    )


def find_computed_table(item):
    """To return the results that already exist"""
    global computed_table, add_find_time, add_hit_time, cont_find_time, cont_hit_time
    if item[0] == "s":
        temp_key = item[1].index_2_key[item[2]]
        the_key = ("s", item[1].weight, item[1].node, temp_key, item[3])
        if computed_table.__contains__(the_key):
            res = computed_table[the_key]
            tdd = TDD(res[1])
            tdd.weight = res[0]
            return tdd
    elif item[0] == "+":
        the_key = (
            "+",
            get_weight(item[1].weight),
            item[1].node,
            get_weight(item[2].weight),
            item[2].node,
        )
        add_find_time += 1
        if computed_table.__contains__(the_key):
            res = computed_table[the_key]
            tdd = TDD(res[1])
            tdd.weight = res[0]
            add_hit_time += 1
            return tdd
        the_key = (
            "+",
            get_weight(item[2].weight),
            item[2].node,
            get_weight(item[1].weight),
            item[1].node,
        )
        if computed_table.__contains__(the_key):
            res = computed_table[the_key]
            tdd = TDD(res[1])
            tdd.weight = res[0]
            add_hit_time += 1
            return tdd
    else:
        the_key = (
            "*",
            get_weight(item[1].weight),
            item[1].node,
            get_weight(item[2].weight),
            item[2].node,
            item[3][0],
            item[3][1],
            item[4],
        )
        cont_find_time += 1
        if computed_table.__contains__(the_key):
            res = computed_table[the_key]
            tdd = TDD(res[1])
            tdd.weight = res[0]
            cont_hit_time += 1
            return tdd
        the_key = (
            "*",
            get_weight(item[2].weight),
            item[2].node,
            get_weight(item[1].weight),
            item[1].node,
            item[3][1],
            item[3][0],
            item[4],
        )
        if computed_table.__contains__(the_key):
            res = computed_table[the_key]
            tdd = TDD(res[1])
            tdd.weight = res[0]
            cont_hit_time += 1
            return tdd
    return None


def insert_2_computed_table(item, res):
    """To insert an item to the computed table"""
    global computed_table, cont_time, find_time, hit_time
    if item[0] == "s":
        temp_key = item[1].index_2_key[item[2]]
        the_key = ("s", item[1].weight, item[1].node, temp_key, item[3])
    elif item[0] == "+":
        the_key = (
            "+",
            get_weight(item[1].weight),
            item[1].node,
            get_weight(item[2].weight),
            item[2].node,
        )
    else:
        the_key = (
            "*",
            get_weight(item[1].weight),
            item[1].node,
            get_weight(item[2].weight),
            item[2].node,
            item[3][0],
            item[3][1],
            item[4],
        )
    computed_table[the_key] = (res.weight, res.node)


def get_index_2_key(var):
    var_sort = copy.copy(var)
    var_sort.sort()
    var_sort.reverse()
    idx_2_key = {-1: -1}
    key_2_idx = {-1: -1}
    n = 0
    for idx in var_sort:
        if not idx.key in idx_2_key:
            idx_2_key[idx.key] = n
            key_2_idx[n] = idx.key
            n += 1
    return idx_2_key, key_2_idx


def get_tdd(U, var=[]):
    idx_2_key, key_2_idx = get_index_2_key(var)
    order = []
    for idx in var:
        order.append(idx_2_key[idx.key])
    tdd = np_2_tdd(U, order)
    tdd.index_2_key = idx_2_key
    tdd.key_2_index = key_2_idx
    tdd.index_set = var

    return tdd


def np_2_tdd(U, order=[], key_width=True):
    # index is the index_set as the axis order of the matrix
    global epi

    if tdd_type == "SymTDD":
        from TDD.SymTDD.BDD import get_bdd
    if tdd_type == "TrDD":
        from TDD.TrDD.BDD import get_bdd
    if tdd_type == "Exp":
        from TDD.Exp.EXP import get_bdd
    U_dim = U.ndim
    U_shape = U.shape
    if sum(U_shape) == U_dim:
        node = Find_Or_Add_Unique_table(-1)
        res = TDD(node)
        for k in range(U_dim):
            U = U[0]

        if isinstance(U, float) and int(round(U / epi)) == 0:
            U = 0
        if isinstance(U, complex):
            r1 = 0
            r2 = 0
            if not int(round(U.real / epi)) == 0:
                r1 = U.real
            if not int(round(U.imag / epi)) == 0:
                r2 = U.imag
            U = complex(r1, r2)

        res.weight = get_bdd(U)

        #         print(res.weight.weight)
        return res

    if not order:
        order = list(range(U_dim))

    if key_width:
        the_width = dict()
        for k in range(max(order) + 1):
            split_pos = order.index(k)
            the_width[k] = U.shape[split_pos]

    x = max(order)
    split_pos = order.index(x)
    order[split_pos] = -1
    split_U = np.split(U, U_shape[split_pos], split_pos)

    while x in order:
        split_pos = order.index(x)
        for k in range(len(split_U)):
            split_U[k] = np.split(split_U[k], U_shape[split_pos], split_pos)[k]
        order[split_pos] = -1

    the_successors = []
    for k in range(U_shape[split_pos]):
        res = np_2_tdd(split_U[k], copy.copy(order), False)
        the_successors.append(res)
    tdd = normalize(x, the_successors)

    if key_width:
        tdd.key_width = the_width

    return tdd


def cont(tdd1, tdd2):
    # 找出哪些要輸出(cont)/保留(out)

    var_cont = [var for var in tdd1.index_set if var in tdd2.index_set]
    var_out1 = [var for var in tdd1.index_set if not var in var_cont]
    var_out2 = [var for var in tdd2.index_set if not var in var_cont]

    var_out = var_out1 + var_out2
    var_out.sort()  # Index 已含有比較大小的函數
    var_out_idx = [var.key for var in var_out]
    var_cont_idx = [var.key for var in var_cont]
    var_cont_idx = [var for var in var_cont_idx if not var in var_out_idx]

    idx_2_key = {-1: -1}
    key_2_idx = {-1: -1}

    n = 0
    for k in range(len(var_out_idx) - 1, -1, -1):
        if not var_out_idx[k] in idx_2_key:
            idx_2_key[var_out_idx[k]] = n
            key_2_idx[n] = var_out_idx[k]
            n += 1

    key_2_new_key = [[], []]  # 找key的對應
    cont_order = [[], []]  # 找cont的順序

    for k in range(len(tdd1.key_2_index) - 1):
        v = tdd1.key_2_index[k]
        if v in idx_2_key:
            key_2_new_key[0].append(idx_2_key[v])
        else:
            key_2_new_key[0].append("c")
        cont_order[0].append(global_index_order[v])

    cont_order[0].append(float("inf"))

    for k in range(len(tdd2.key_2_index) - 1):
        v = tdd2.key_2_index[k]
        if v in idx_2_key:
            key_2_new_key[1].append(idx_2_key[v])
        else:
            key_2_new_key[1].append("c")
        cont_order[1].append(global_index_order[v])
    cont_order[1].append(float("inf"))

    tdd = contract(tdd1, tdd2, key_2_new_key, cont_order, len(set(var_cont_idx)))
    tdd.index_set = var_out
    tdd.index_2_key = idx_2_key
    tdd.key_2_index = key_2_idx
    key_width = dict()
    for k1 in range(len(key_2_new_key[0])):
        if not key_2_new_key[0][k1] == "c" and not key_2_new_key[0][k1] == -1:
            key_width[key_2_new_key[0][k1]] = tdd1.key_width[k1]
    for k2 in range(len(key_2_new_key[1])):
        if not key_2_new_key[1][k2] == "c" and not key_2_new_key[1][k2] == -1:
            key_width[key_2_new_key[1][k2]] = tdd2.key_width[k2]

    tdd.key_width = key_width
    #     print(tdd1.key_width,tdd2.key_width,tdd.key_width)
    return tdd


def mul_weight(w1, w2):
    res = w1 * w2
    return res


def add_weight(w1, w2):
    #     t1=time.time()
    res = w1 + w2
    #     t=time.time()-t1
    #     if t>0.001:
    #         print('add',t)
    #         print(w1)
    #         print(w2)
    return res


def contract(tdd1, tdd2, key_2_new_key, cont_order, cont_num):
    """The contraction of two TDDs, var_cont is in the form [[4,1],[3,2]]"""
    global S_one, S_zero

    k1 = tdd1.node.key
    k2 = tdd2.node.key
    w1 = tdd1.weight
    w2 = tdd2.weight

    if k1 == -1 and k2 == -1:
        if w1 == S_zero:
            tdd = TDD(tdd1.node)
            tdd.weight = S_zero
            return tdd
        if w2 == S_zero:
            tdd = TDD(tdd1.node)
            tdd.weight = S_zero
            return tdd
        tdd = TDD(tdd1.node)
        tdd.weight = mul_weight(w1, w2)
        if cont_num > 0:
            tdd.weight.weight *= 2**cont_num
        return tdd
    # 注意，我们在这里允许直接乘上一个权重，而不是把它push到底层去，可能会导致表示不唯一，若要进行等价性检验，需要在最后renormalise一下
    if k1 == -1:
        if w1 == S_zero:
            tdd = TDD(tdd1.node)
            tdd.weight = S_zero
            return tdd
        if cont_num == 0 and key_2_new_key[1][k2] == k2:
            tdd = TDD(tdd2.node)
            tdd.weight = mul_weight(w1, w2)
            if tdd.weight == S_zero:
                tdd.node = Find_Or_Add_Unique_table(-1)
            return tdd

    if k2 == -1:
        if w2 == S_zero:
            tdd = TDD(tdd2.node)
            tdd.weight = S_zero
            return tdd
        if cont_num == 0 and key_2_new_key[0][k1] == k1:
            tdd = TDD(tdd1.node)
            tdd.weight = mul_weight(w1, w2)
            if tdd.weight == S_zero:
                tdd.node = Find_Or_Add_Unique_table(-1)
            return tdd

    tdd1.weight = S_one
    tdd2.weight = S_one
    #     print('ccc:',S_one.node)
    temp_key_2_new_key = []
    temp_key_2_new_key.append(tuple([k for k in key_2_new_key[0][: (k1 + 1)]]))
    temp_key_2_new_key.append(tuple([k for k in key_2_new_key[1][: (k2 + 1)]]))

    tdd = find_computed_table(["*", tdd1, tdd2, temp_key_2_new_key, cont_num])
    if tdd:
        #         tdd.weight=simplify(tdd.weight*w1*w2)
        tdd.weight = mul_weight(tdd.weight, mul_weight(w1, w2))
        tdd1.weight = w1
        tdd2.weight = w2
        if tdd.weight == S_zero:
            tdd.node = Find_Or_Add_Unique_table(-1)
        return tdd

    if cont_order[0][k1] < cont_order[1][k2]:
        the_key = key_2_new_key[0][k1]
        if the_key != "c":
            the_successors = []
            for k in range(tdd1.node.succ_num):
                res = contract(
                    Slicing(tdd1, k1, k), tdd2, key_2_new_key, cont_order, cont_num
                )
                the_successors.append(res)
            tdd = normalize(the_key, the_successors)
            insert_2_computed_table(
                ["*", tdd1, tdd2, temp_key_2_new_key, cont_num], tdd
            )
            tdd.weight = mul_weight(tdd.weight, mul_weight(w1, w2))
        else:
            tdd = TDD(Find_Or_Add_Unique_table(-1))
            tdd.weight = S_zero
            for k in range(tdd1.node.succ_num):
                res = contract(
                    Slicing(tdd1, k1, k), tdd2, key_2_new_key, cont_order, cont_num - 1
                )
                tdd = add(tdd, res)
            insert_2_computed_table(
                ["*", tdd1, tdd2, temp_key_2_new_key, cont_num], tdd
            )
            tdd.weight = mul_weight(tdd.weight, mul_weight(w1, w2))
    elif cont_order[0][k1] == cont_order[1][k2]:
        the_key = key_2_new_key[0][k1]
        if the_key != "c":
            the_successors = []
            for k in range(tdd1.node.succ_num):
                res = contract(
                    Slicing(tdd1, k1, k),
                    Slicing(tdd2, k2, k),
                    key_2_new_key,
                    cont_order,
                    cont_num,
                )
                the_successors.append(res)
            tdd = normalize(the_key, the_successors)
            insert_2_computed_table(
                ["*", tdd1, tdd2, temp_key_2_new_key, cont_num], tdd
            )
            tdd.weight = mul_weight(tdd.weight, mul_weight(w1, w2))
        else:
            tdd = TDD(Find_Or_Add_Unique_table(-1))
            tdd.weight = S_zero
            for k in range(tdd1.node.succ_num):
                res = contract(
                    Slicing(tdd1, k1, k),
                    Slicing(tdd2, k2, k),
                    key_2_new_key,
                    cont_order,
                    cont_num - 1,
                )
                tdd = add(tdd, res)
            insert_2_computed_table(
                ["*", tdd1, tdd2, temp_key_2_new_key, cont_num], tdd
            )
            tdd.weight = mul_weight(tdd.weight, mul_weight(w1, w2))
    else:
        the_key = key_2_new_key[1][k2]
        if the_key != "c":
            the_successors = []
            for k in range(tdd2.node.succ_num):
                res = contract(
                    tdd1, Slicing(tdd2, k2, k), key_2_new_key, cont_order, cont_num
                )
                the_successors.append(res)
            tdd = normalize(the_key, the_successors)
            insert_2_computed_table(
                ["*", tdd1, tdd2, temp_key_2_new_key, cont_num], tdd
            )
            tdd.weight = mul_weight(tdd.weight, mul_weight(w1, w2))
        else:
            tdd = TDD(Find_Or_Add_Unique_table(-1))
            tdd.weight = S_zero
            for k in range(tdd2.node.succ_num):
                res = contract(
                    tdd1, Slicing(tdd2, k2, k), key_2_new_key, cont_order, cont_num - 1
                )
                tdd = add(tdd, res)
            insert_2_computed_table(
                ["*", tdd1, tdd2, temp_key_2_new_key, cont_num], tdd
            )
            tdd.weight = mul_weight(tdd.weight, mul_weight(w1, w2))
    tdd1.weight = w1
    tdd2.weight = w2
    #     tdd.weight=simplify(tdd.weight)
    if tdd.weight == S_zero:
        tdd.node = Find_Or_Add_Unique_table(-1)
    return tdd


def Slicing(tdd, x, c):
    """Slice a TDD with respect to x=c"""

    k = tdd.node.key

    if k == -1:
        return tdd.self_copy()

    if k < x:
        return tdd.self_copy()

    if k == x:
        res = TDD(tdd.node.successor[c])
        res.weight = tdd.node.out_weight[c]
        return res
    else:
        print("Not supported yet!!!")


def Slicing2(tdd, x, c):
    """Slice a TDD with respect to x=c"""

    k = tdd.node.key

    if k == -1:
        return tdd.self_copy()

    if k < x:
        return tdd.self_copy()

    if k == x:
        res = TDD(tdd.node.successor[c])
        #         res.weight=simplify(tdd.node.out_weight[c]*tdd.weight)
        res.weight = mul_weight(tdd.weight, tdd.node.out_weight[c])
        return res
    else:
        print("Not supported yet!!!")


def add(tdd1, tdd2):
    """The apply function of two TDDs. Mostly, it is used to do addition here."""
    global global_index_order, S_one, S_zero

    k1 = tdd1.node.key
    k2 = tdd2.node.key

    if tdd1.weight == S_zero:
        return tdd2.self_copy()

    if tdd2.weight == S_zero:
        return tdd1.self_copy()

    if tdd1.node == tdd2.node:
        weig = add_weight(tdd1.weight, tdd2.weight)
        if weig == S_zero:
            term = Find_Or_Add_Unique_table(-1)
            res = TDD(term)
            res.weight = S_zero
            return res
        else:
            res = TDD(tdd1.node)
            res.weight = weig
            return res

    if find_computed_table(["+", tdd1, tdd2]):
        return find_computed_table(["+", tdd1, tdd2])
    the_successors = []
    if k1 > k2:
        x = k1
        for k in range(tdd1.node.succ_num):
            res = add(Slicing2(tdd1, x, k), tdd2)
            the_successors.append(res)
    elif k1 == k2:
        x = k1
        for k in range(tdd1.node.succ_num):
            res = add(Slicing2(tdd1, x, k), Slicing2(tdd2, x, k))
            the_successors.append(res)
    else:
        x = k2
        for k in range(tdd2.node.succ_num):
            res = add(tdd1, Slicing2(tdd2, x, k))
            the_successors.append(res)

    res = normalize(x, the_successors)
    insert_2_computed_table(["+", tdd1, tdd2], res)
    return res


def conjugate2(tdd):
    """ "change the position of x and y, x<y in this tdd"""
    v = tdd.node.key
    if v == -1:
        res = tdd.self_copy()
        res.weight = tdd.weight.conj()
        return res
    low = conjugate2(Slicing2(tdd, v, 0))
    high = conjugate2(Slicing2(tdd, v, 1))
    res = normalize(v, [low, high])
    res.index_set = [copy.copy(k) for k in tdd.index_set]
    res.key_2_index = copy.copy(tdd.key_2_index)
    res.index_2_key = copy.copy(tdd.index_2_key)
    res.key_width = copy.copy(tdd.key_width)
    return res


def global_norm(tdd):
    node = tdd.node
    weight = tdd.weight

    if node.key == -1:
        return tdd

    h = weight
    tdd_list = []
    for i in range(len(node.out_weight)):
        tdd_list.append(TDD(node.successor[i]))
        tdd_list[i].weight = h * node.out_weight[i]
        tdd_list[i] = global_norm(tdd_list[i])

    new_tdd = normalize(node.key, tdd_list)
    return new_tdd
