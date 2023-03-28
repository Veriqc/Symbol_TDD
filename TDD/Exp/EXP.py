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
from functools import lru_cache, wraps, update_wrapper

"""Define global variables"""

unique_table = dict()

global_index_order = dict()

var_num = 0

epi=1e-5


class BDD:
    def __init__(self,data=dict()):
        """An expression of c exp(\sum{b_i \theta_i})"""
        self.data = data

    @property
    def node(self):
        """Convert to hashable tuple"""

        return tuple(sorted(self.data.items(), key=lambda x: x[0]))

    @property
    def weight(self):
        """Reserved to support current cache table key construction."""
        return 0
    
    @property
    def is_constant(self):
        return not all( next(iter(self.data.keys())) ) if len(self.data) == 1 else False
    
    @property
    def is_zero(self):
        return self == BDD()

    def __repr__(self) -> str:
        return str(sorted(self.data.items(), key=lambda x: x[0]))
    
    def __hash__(self):
        return hash(self.__repr__())

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
        print('EXP 75', item)
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
    if theta < 0:  # 將小於0的角度+2pi
        theta += 2*np.pi
        # print(theta)
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

def get_bdd(function):
    global global_index_order

    if isinstance(function, BDD):
        return function
    
    elif not isinstance(function, sp.Basic) or len(function.args) == 0:
        return get_const_state(function)
    
    elif len(function.args) == 1:
        # print(function.args[0])
        def get_item(expr):
            dict1 = dict()
            dict2 = dict()
            for item in expr.free_symbols:
                dict1[item] = expr.coeff(item)/1j
                dict2[item] = 0
            return dict1, float(expr.subs(dict2)/1j)
        
        dict1, theta = get_item(function.args[0])
        data = [0]*(var_num)

        for key in dict1.keys():
            data[global_index_order[key]] = dict1[key]

        return BDD({tuple(data):(1, theta)})
    
    if isinstance(function, sp.core.add.Add):
        tdd = BDD()
        for add_term in function.args:
            temp_tdd = get_bdd(add_term)
            tdd = add(tdd, temp_tdd)
    elif isinstance(function, sp.core.mul.Mul):
        tdd = get_one_state()
        for mul_term in function.args:
            temp_tdd = get_bdd(mul_term)
            tdd = mul(tdd, temp_tdd)
    elif isinstance(function, sp.core.power.Pow):
        tdd = get_one_state()
        pow = function.args[1]
        for mul_times in range(pow):
            temp_tdd = get_bdd(function.args[0])
            tdd = mul(tdd, temp_tdd)
    return tdd

def add_constant(const1:tuple, const2:tuple):
    value = complex(const1[0]*exp(1j*const1[1])+const2[0]*exp(1j*const2[1]))
    r, theta = complex_to_polar(value)
    
    return (np.round(r,int(-np.log10(epi))), np.round(theta,int(-np.log10(epi))))

def hashed_args(func):
    @wraps(func)
    def wrapper(*args):
        # print(*tuple(sorted(args)))
        return func(*tuple(sorted(args, key= hash)))
    
    update_wrapper(wrapper, func, ('cache_info', 'cache_clear','cache_parameters')) #繼承lru_cache的屬性
    return wrapper

@hashed_args #為了達成交換率，先排序輸入的兩個bdd
@lru_cache(maxsize=None)
def add(bdd1, bdd2):
    """The apply function of two TDDs. Mostly, it is used to do addition here."""
    
    keys_intersect = bdd1.data.keys() & bdd2.data.keys()

    data = copy.copy(bdd1.data)
    data.update(bdd2.data)
    # print('EXP 175' , data)
    for term in keys_intersect:
        value = add_constant(bdd1.data[term], bdd2.data[term])
        if not math.isclose(value[0], 0 , abs_tol=epi):
            data[term] = value 
        else:
            data.pop(term)
    # print('EXP 182',data)
    return BDD(data)

@hashed_args #為了達成交換率，先排序輸入的兩個bdd
@lru_cache(maxsize=None)
def mul(bdd1, bdd2):
    data=dict()
    for term1 in bdd1.data:
        for term2 in bdd2.data:
            # t=tuple( term1[i]+term2[i] for i in range(len(term1)))

            t=tuple( np.round(float(term1[i]+term2[i]),int(-np.log10(epi))) for i in range(len(term1)))

            new_r = np.round(bdd1.data[term1][0] * bdd2.data[term2][0],int(-np.log10(epi)))
            new_theta = np.round(bdd1.data[term1][1] + bdd2.data[term2][1],int(-np.log10(epi)))
            if new_theta > 2 * np.pi: #將>2pi的角度-2pi
                new_theta -= 2 * np.pi 
            if math.isclose(new_theta, 2*np.pi, abs_tol=epi): #將=2pi角度改為0.0
                new_theta = 0.0

            if t in data:
                new_r, new_theta = add_constant(data[t],(new_r,new_theta))

            if not math.isclose(new_r, 0, abs_tol=epi): 
                data[t] = (new_r, new_theta) 
            else:
                data.pop(t, None) 
    return BDD(data)

def normalize_2_fun(bdd1,bdd2):

    return [get_one_state(),bdd1,bdd2]       
    
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
