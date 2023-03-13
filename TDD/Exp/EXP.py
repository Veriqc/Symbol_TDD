import numpy as np
import copy
import time
import random
from graphviz import Digraph
from IPython.display import Image
from sympy import *
from sympy.parsing.sympy_parser import parse_expr
import sympy as sp
import math


"""Define global variables"""

unique_table = dict()

global_index_order = dict()

var_num = 0

epi=0.00000001


class BDD:
    def __init__(self,data=dict()):
        """An expression of c exp(\sum{b_i \theta_i})"""
        self.data = data

    @property
    def node(self):
        """Convert to hashable tuple"""
        # print(tuple(self.data.items()))
        return tuple(self.data.items())

    @property
    def weight(self):
        """Reserved to support current cache table key construction."""
        return 0
    
    @property
    def is_constant(self):
        return all( next(iter(self.data.keys())) ) if len(self.data) == 1 else False
    
    @property
    def is_zero(self):
        return self == BDD()

    def __repr__(self) -> str:
        return str(self.data)

    def __eq__(self,other):
        return self.data == other.data
        
    def __add__(self, g):        
        return add(self, g)
    
    def __mul__(self, g):
        return mul(self, g)
    
    # def __div__(self, g):
    #     return normalize_2_fun(g, self)
    
    def __sub__(self, g):
        return self.__add__(g.__mul__(-1))
    
    def self_normalize(self, g):
        # self/g
        if g.is_zero:
            return self, g, get_one_state()

        key1=min(g.data.keys())

        item={key1:g.data[key1]}

        def div_item (data,item ):
            dict1={}
            key=list(item.keys())[0]
            for K in data.keys():
                new_key=tuple([K[i] - key[i] for i in range(len(key))])
                value = (data[K][0]/item[key][0],data[K][1]-item[key][1] ) 
                dict1[new_key]=value
            return dict1

        r = div_item(self.data, item)
        l = div_item(g.data, item)

        return BDD(item), BDD(l) , BDD(r)
    
        
def Ini_BDD(index_order=[]):
    """To initialize the unique_table, computed_table and set up a global index order"""
    global global_index_order, var_num
    
    set_index_order(index_order)

    var_num = len(index_order)

def complex_to_polar(c):
    r = np.abs(c)
    theta = np.angle(c)
    return r, theta

def get_one_state():
    #key 是符號項係數，value 是 [r,theta]
    data = {tuple([0]*var_num):(1, 0)}
    return BDD(data)

def get_const_state(value):
    #key 是符號項係數，value 是 [r,theta]
    value = complex(value)
    r, theta = complex_to_polar(value)
    key = [0]*(var_num)
    # data = {tuple(key):(r, theta)}
    return BDD({} if math.isclose(r, 0, abs_tol=epi) else {tuple(key):(r, theta)})

def set_index_order(var_order):
    global global_index_order
    global_index_order=dict()
    if isinstance(var_order,list):
        for k in range(len(var_order)):
            global_index_order[var_order[k]]=k
    if isinstance(var_order,dict):
        global_index_order = copy.copy(var_order)
    global_index_order[1] = var_num
    
def get_index_order():
    global global_index_order
    return copy.copy(global_index_order)
    
def get_int_key(weight):
    """To transform a complex number to a tuple with int values"""
    global epi
    return ( int(round(weight.real/epi)), int(round(weight.imag/epi)) )

def get_bdd(f):
    global global_index_order

    if isinstance(f, BDD):
        return f
    
    elif not isinstance(f, sp.Basic) or len(f.args) == 0:
        return get_const_state(f)
    elif len(f.args) == 1:
        # print(f.args[0])
        def get_item(expr):
            dict1 = dict()
            dict2 = dict()
            for item in expr.free_symbols:
                dict1[item] = expr.coeff(item)/1j
                dict2[item] = 0
            return dict1, float(expr.subs(dict2)/1j)
        
        dict1, theta = get_item(f.args[0])
        data = [0]*(var_num)

        for key in dict1.keys():
            data[global_index_order[key]] = dict1[key]

        return BDD({tuple(data):(1, theta)})
    
    if isinstance(f, sp.core.add.Add):
        tdd = BDD()
        for add_term in f.args:
            temp_tdd = get_bdd(add_term)
            tdd = add(tdd, temp_tdd)
    elif isinstance(f, sp.core.mul.Mul):
        tdd = get_one_state()
        for mul_term in f.args:
            temp_tdd = get_bdd(mul_term)
            tdd = mul(tdd, temp_tdd)
    elif isinstance(f, sp.core.power.Pow):
        tdd = get_one_state()
        pow = f.args[1]
        for mul_times in range(pow):
            temp_tdd = get_bdd(f.args[0])
            tdd = mul(tdd, temp_tdd)
    return tdd

def add_constant(const1:tuple, const2:tuple):
    value = complex(const1[0]*exp(1j*const1[1])+const2[0]*exp(1j*const2[1]))
    r, theta = complex_to_polar(value)
    
    return (np.round(r,int(-np.log10(epi))), np.round(theta,int(-np.log10(epi))))

def add(tdd1, tdd2):
    """The apply function of two TDDs. Mostly, it is used to do addition here."""
    
    ### Necessary?
    # if len(tdd2.data) > len(tdd1.data):
    #     return add(tdd2,tdd1)
    
    keys_intersect = tdd1.data.keys() & tdd2.data.keys()

    data = copy.copy(tdd1.data)
    data.update(tdd2.data)
    # print('EXP 175' , data)
    for term in keys_intersect:
        value = add_constant(tdd1.data[term], tdd2.data[term])
        if not math.isclose(value[0], 0 , abs_tol=epi):
            data[term] = value 
        else:
            data.pop(term)
    # print('EXP 182',data)
    return BDD(data)

def mul(tdd1, tdd2):
    # print('191',tdd1.data)
    # print('192',tdd2.data)
    """The contraction of two TDDs, var_cont is in the form [[4,1],[3,2]]"""
    # if not isinstance(tdd1,BDD): 
    #     data=tdd2.data.copy()
    #     for k in data:
    #         data[k]*=tdd1
    #     return BDD(data)
    # if not isinstance(tdd2,BDD):
    #     data=tdd1.data.copy()
    #     for k in data:
    #         data[k]*=tdd2
    #     return BDD(data)
    
    data=dict()
    
    for term1 in tdd1.data:
        for term2 in tdd2.data:
            t=tuple(term1[i]+term2[i] for i in range(len(term1)))
            new_r = tdd1.data[term1][0]*tdd2.data[term2][0]
            new_theta = tdd1.data[term1][1]+tdd2.data[term2][1]

            if t in data:
                new_r, new_theta = add_constant(data[t],(new_r,new_theta))

            if not math.isclose(new_r, 0, abs_tol=epi): 
                data[t] = (new_r, new_theta) 
            else:
                data.pop(t) 
    # print('EXP 232',data)
    return BDD(data)
        




def normalize_2_fun(tdd1,tdd2):

    return [get_one_state(),tdd1,tdd2]      
    
    
    
    
def get_expr(node):
    if node.key==-1:
        return 1
    if node.expr:
        return node.expr
    x=node.key
    
    expr=parse_expr(x)
    
    res=0
    for succ in node.successor:
        res+=succ[1]*expr**succ[0]*get_expr(succ[2])

    node.expr=res
    return res

def get_value_node(node,val,the_key=None):
    if not the_key:
        the_key=tuple([tuple([k,val[k]]) for k in val])
    if the_key in node.value:
        return node.value[the_key]
    
    if node.key==-1:
        return 1
    
    res=0
    for succ in node.successor:
        expr=parse_expr(node.key)
        v=complex(expr.subs(val))
        res+=succ[1]*v**succ[0]*get_value_node(succ[2],val,the_key)

    node.value[the_key]=res
    return res

def get_value_node2(node,val,the_key=None):
    if not the_key:
        the_key=tuple([tuple([k,val[k]]) for k in val])
    if the_key in node.value:
        return node.value[the_key]
    
    if node.key==-1:
        return 1
    
    res=0

    v=val[node.key]
    for succ in node.successor:
        res+=succ[1]*(v**succ[0])*get_value_node2(succ[2],val,the_key)
    node.value[the_key]=res
    return res

def conjugate(tdd):
    """"change the position of x and y, x<y in this tdd
    """
    v=tdd.node.key
    if v==-1:
        res=tdd.self_copy()
        res.weight=tdd.weight.conjugate()      
        return res
    the_successor=[]
    for k in tdd.node.successor:
        temp_res=conjugate(BDD(k[2]))
        the_successor.append([k[0],k[1].conjugate()*temp_res.weight,temp_res.node])
    res=normalize(v,the_successor)
    return res
