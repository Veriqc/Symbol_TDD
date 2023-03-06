import numpy as np
import copy
import time
import random
from graphviz import Digraph
from IPython.display import Image
from sympy import *
from sympy.parsing.sympy_parser import parse_expr
import sympy as sy


"""Define global variables"""

unique_table = dict()

global_index_order = dict()

var_num = 0

epi=0.00000001


class BDD:
    def __init__(self,data=None):
        """An expression of c exp(\sum{b_i \theta_i})"""
        if data:
            self.data=data
        else:
            self.data=dict()
    @property
    def node(self):
        # print(self.data,type(self))
        from functools import reduce
        return reduce(lambda a,b:a+b , self.data.keys())
    @property
    def weight(self):
        return 0

    def __eq__(self,other):
        if self.data==other.data:
            return True
        else:
            return False
        
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
    
        
def Ini_BDD(index_order=[],Exp_parameters_dict={}):
    """To initialize the unique_table,computed_table and set up a global index order"""
    global global_node_idx, var_num, parameters_dict, parameter_num
    
    
    set_index_order(index_order)
    var_num = len(index_order)
    parameters_dict=Exp_parameters_dict
    parameter_num=len(parameters_dict)

def get_one_state():
    data={tuple([0]*(parameter_num+1)):1}
    tdd = BDD(data)
    return tdd

def get_zero_state():
    data={tuple([0]*(parameter_num+1)):0}
    tdd = BDD(data)
    return tdd

def get_const_state(w):
    r=np.abs(w)
    theta = np.angle(w)
    key = [0]*(parameter_num+1)
    key[-1] = theta
    data = {tuple(key):r}
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
    global_index_order[-1] = float('inf')
    
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
    print(f)

    def extract_expr(param):
        if len(param.parameters)==0:
            return None, None, float(param)
        p=tuple(param.parameters)[0]
        const=param.bind({p:0})
        const=const.sympify().simplify()
        coeff=param.bind({p:1}).sympify()-const
        coeff=coeff.simplify()
        return float(coeff), p, float(const)
    # expr= [extract_expr(item) for item in g[0].params]

    def construct_exp_expr(coeff,p,const):
        r=np.abs(const)
        theta = np.angle(const)
        temp=[0]*parameter_num+[theta]
        if p:
            pos=parameters_dict[p]
            temp[pos]=coeff
        from TDD.Exp.EXP import BDD
        return BDD({tuple(temp):r})
    
    def angle_mul(theta, scale):
        return [theta[0]*scale,theta[1],theta[2]*scale]

    return BDD(f)


def add(tdd1,tdd2):
    """The apply function of two TDDs. Mostly, it is used to do addition here."""
    
    if len(tdd2.data)>len(tdd1.data):
        return add(tdd2,tdd1)
    
    data=copy.copy(tdd1.data)
    
    for term in tdd2.data:
        if term in data:
            data[term]+=tdd2.data[term]
        else:
            data[term]=tdd2.data[term]
    
    return BDD(data)



def mul(tdd1,tdd2):
    """The contraction of two TDDs, var_cont is in the form [[4,1],[3,2]]"""
    if not isinstance(tdd1,BDD): 
        data=tdd2.data.copy()
        for k in data:
            data[k]*=tdd1
        return BDD(data)
    if not isinstance(tdd2,BDD):
        data=tdd1.data.copy()
        for k in data:
            data[k]*=tdd2
        return BDD(data)
    
    data=dict()
    
    for term1 in tdd1.data:
        for term2 in tdd2.data:
            t=[term1[i]+term2[i] for i in range(len(term1))]
            t=tuple(t)
            if t in data:
                data[t]+=tdd1.data[term1]*tdd2.data[term2]
            else:
                data[t]=tdd1.data[term1]*tdd2.data[term2]

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
