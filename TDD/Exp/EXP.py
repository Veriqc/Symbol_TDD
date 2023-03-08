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
    def __init__(self,data=None):
        """An expression of c exp(\sum{b_i \theta_i})"""
        self.data=data if data else dict()
    @property
    def node(self):
        # print(self.data,type(self))
        # from functools import reduce
        l=[]
        for k in self.data.keys():
            l+=list(k)
            l+=self.data[k]
        return tuple(l)
        
        # return reduce(lambda a,b:a+b , self.data.keys())

    @property
    def weight(self):
        return 0
    @property
    def is_constant(self):
        keys = list(self.data.keys())
        return len(keys)==1 and all(keys[0][i]==0 for i in range (len(keys[0])-1))
    @property
    def is_zero(self):
        return self==BDD()

    def __repr__(self) -> str:
        return str(self.data)

    def __eq__(self,other):
        return self.data==other.data
        
    def __add__(self, g):        
        return add(self,g)
    
    def __mul__(self, g):
        return mul(self,g)
    
    def __div__(self, g):
        return normalize_2_fun(g,self)
    def __sub__(self, g):
        return self.__add__(g.__mul__(-1))
    
    def self_normalize(self, g):
        return get_one_state(),g , self 
    
        
def Ini_BDD(index_order=[]):
    """To initialize the unique_table,computed_table and set up a global index order"""
    global global_index_order, var_num
    
    # print('EXP 63', index_order)
    set_index_order(index_order)
    # print('EXP 65', global_index_order)
    var_num = len(index_order)


def get_one_state():
    #key 是符號項係數，value 是 [r,theta]
    data={tuple([0]*var_num):[1,0]}
    tdd = BDD(data)
    return tdd


def get_const_state(w):
    #key 是符號項係數，value 是 [r,theta]
    w=complex(w)
    # print('EXP 79',w)
    r=np.abs(w)
    theta = np.angle(w)
    key = [0]*(var_num)
    data = {tuple(key):[r,theta]}
    tdd = BDD(data)
    return tdd

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
    return (int(round(weight.real/epi)) ,int(round(weight.imag/epi)))


def get_bdd(f):
    global global_index_order 
    if isinstance(f,int) or isinstance(f,float) or isinstance(f,complex):
        tdd=get_const_state(f)
        return tdd
    if isinstance(f,BDD):
        return f
    
    if len(f.args)==0:
        tdd=get_const_state(f)
        return tdd
    
    if len(f.args)==1:
        def get_item(f):
            expr=f.args[0]
            dict1=dict()
            dict2=dict()
            for item in expr.free_symbols:
                dict1[item]=expr.coeff(item)
                dict2[item]=0
            # dict1[1]=expr.subs(dict2)
            return dict1, complex(expr.subs(dict2))
        dict1,c =get_item(f)
        data=[0]*(var_num)

        for key in dict1.keys():
            data[global_index_order[key]]=dict1[key]
        dict2={tuple(data):[1,c]}
        tdd = BDD(dict2)
        return tdd
    if isinstance(f,sp.core.add.Add):
        tdd=BDD()
        for add_term in f.args:
            temp_tdd=get_bdd(add_term)
            tdd=add(tdd,temp_tdd)
    elif isinstance(f,sp.core.mul.Mul):
        tdd=get_one_state()
        for mul_term in f.args:
            temp_tdd=get_bdd(mul_term)
            tdd=mul(tdd,temp_tdd)
    elif isinstance(f,sp.core.power.Pow):
        tdd=get_one_state()
        pow= f.args[1]
        for mul_times in range(pow):
            temp_tdd=get_bdd(f.args[0])
            tdd=mul(tdd,temp_tdd)
    return tdd

def add_constant(const1:list,const2:list):
    value=complex(const1[0]*exp(1j*const1[1])+const2[0]*exp(1j*const2[1]))
    r=np.abs(value)
    theta = np.angle(value)

    return [r,theta]
def add(tdd1,tdd2):
    """The apply function of two TDDs. Mostly, it is used to do addition here."""
    
    if len(tdd2.data)>len(tdd1.data):
        return add(tdd2,tdd1)
    
    data=copy.copy(tdd1.data)
    
    for term in tdd2.data:
        if term in data:
            data[term]=add_constant(tdd1.data[term],tdd2.data[term])
            if math.isclose(data[term][0],epi):
                data.pop(term,None)
        else:
            data[term]=tdd2.data[term]
    
    return BDD(data)



def mul(tdd1,tdd2):
    # print('191',tdd1.data)
    # print('192',tdd2.data)
    """The contraction of two TDDs, var_cont is in the form [[4,1],[3,2]]"""
    if not isinstance(tdd1,BDD): 
        data=tdd2.data.copy()
        for k in data:
            data[k]*=tdd1
        # print('EXP 196',data)
        return BDD(data)
    if not isinstance(tdd2,BDD):
        data=tdd1.data.copy()
        for k in data:
            data[k]*=tdd2
        # print('EXP 202',data)
        return BDD(data)
    
    data=dict()
    
    for term1 in tdd1.data:
        not_pop1=True
        if math.isclose(tdd1.data[term1][0],0): #pop the zero terms
            tdd1.data.copy().pop(term1,None)
            not_pop1=False
        for term2 in tdd2.data:
            not_pop2=True
            if math.isclose(tdd2.data[term2][0],0): #pop the zero terms
                tdd2.data.copy().pop(term2,None)
                not_pop2=False
            t=[term1[i]+term2[i] for i in range(len(term1))]
            t=tuple(t)
            # print('EXP 207', t)
            if not_pop1 and not_pop2:
                if t in data:
                    data[t][0]*=tdd1.data[term1][0]*tdd2.data[term2][0]
                    data[t][1]+=tdd1.data[term1][1]+tdd2.data[term2][1]
                else:
                    data[t]=[tdd1.data[term1][0]*tdd2.data[term2][0],tdd1.data[term1][1]+tdd2.data[term2][1]]

    print('EXP 214',data)
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
