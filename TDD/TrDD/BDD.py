import numpy as np
import copy
import math
from graphviz import Digraph
from IPython.display import Image
from sympy import *
from sympy.parsing.sympy_parser import parse_expr
import sympy as sp
from functools import lru_cache

"""Define global variables"""
computed_table = dict()
unique_table = dict()
global_index_order = dict()
global_node_idx=0
add_find_time=0
add_hit_time=0
cont_find_time=0
cont_hit_time=0
epi=1e-10

    
class Node:
    """To define the node of TDD"""
    def __init__(self,key):
        self.idx = 0
        self.key = key
        self.successor=[]#the element should be in this form [degree,weight,node]
        self.index_degree=dict()
        self.expr = None
        self.value=dict()
        self.succ_num=2

    def get_expr(self):
        if self.key==-1:
            return 1
        if self.expr:
            return self.expr
        
        def get_numbers(s):
            import re
            # 使用正則表達式找到所有匹配"\d+"的子串，即連續的一個或多個數字
            numbers = re.findall("\d+", s)
            return int(numbers[0])
        param_expr=self.key
        sym_str=param_expr.replace('sin(','').replace('cos(','').replace('[','').replace('])','').replace(str(get_numbers(param_expr)),'')

        sym_base = IndexedBase(sym_str)

        sp_expr=sympify(str(param_expr), locals={sym_str: sym_base})
        

        res=0
        for succ in self.successor:
            res+=nsimplify(succ[1]*sp_expr**succ[0],tolerance=1e-3)*succ[2].get_expr()

        self.expr=res
        return res


class BDD:
    def __init__(self,node):
        """BDD"""
        self._weight=1

        self.key_2_index=dict()

        if isinstance(node,Node):
            self.node=node
        else:
            self.node=Node(node)

    @property
    def weight(self):
        return self._weight
    @weight.setter
    def weight(self, value):
        value=[value.real,value.imag]

        for i in range(2):
            if math.isclose(value[i] , int(value[i]), rel_tol = epi):
                value[i] = int(value[i])
            elif math.isclose(value[i] , int(value[i]+1), rel_tol = epi):
                value[i] = int(value[i]+1)
            elif math.isclose(value[i]+1 , int(value[i]+1), rel_tol = epi):
                value[i] = int(value[i]+1)-1

        self._weight=value[0]+value[1]*1j

    def node_number(self):
        node_set=set()
        node_set=get_node_set(self.node,node_set)
        return len(node_set)
    
    def self_copy(self):
        temp = BDD(self.node)
        temp.weight = self.weight
        return temp
    
    def show(self,real_label=True):
        edge=[]              
        dot=Digraph(name='reduced_tree')
        dot=layout(self.node,dot,edge)
        dot.node('o','',shape='none')
        # dot.edge('-0',str(self.node.idx),color="blue",label=str(complex(round(self.weight.real,2),round(self.weight.imag,2))))
        dot.edge('o',str(self.node.idx),color="blue",label=str(nsimplify(self.weight,rational=False,tolerance=1e-3)))
        dot.format = 'png'
        return Image(dot.render('output'))
    def conj(self):
        w=self.weight
        if w==0:
            return get_zero_state()
        self.weight=1
        res=conjugate(self)
        res.weight*=w.conjugate()
        self.weight=w
        return res
        
    def __eq__(self,other):
        return self.node==other.node and math.isclose(self.weight.real, other.weight.real, rel_tol=epi) and math.isclose(self.weight.imag, other.weight.imag, rel_tol=epi)
        
    def __add__(self, g):        
        return add(self,g)
    
    def __mul__(self, g):
        return mul(self,g)
    
    def self_normalize(self, g):
        return normalize_2_fun(g,self)
    
    # @lru_cache(maxsize=None)
    def expr(self):

        value=[self.weight.real,self.weight.imag]
        for i in range(2):
            if math.isclose(value[i] , int(value[i]), rel_tol = epi):
                value[i] = int(value[i])
            elif math.isclose(value[i] , int(value[i]+1), rel_tol = epi):
                value[i] = int(value[i]+1)
            elif math.isclose(value[i]+1 , int(value[i]+1), rel_tol = epi):
                value[i] = int(value[i]+1)-1
        value=value[0]+value[1]*I
        if self.node.key==-1:
            return value
        
        res=self.node.get_expr()
        
        return  value*res
    
    def get_value(self,val):
        res=get_value_node(self.node,val)
        return self.weight*res
    
    def get_value2(self,val):
        res=get_value_node2(self.node,val)
        return self.weight*res    
    
    def __hash__(self):
        return hash(self.__repr__())
    
    def __repr__(self):
        return str(self.expr())
        
    
def layout(node,dot=Digraph(),succ=[]):
    col=['red','blue','black','green']
    dot.node(str(node.idx), str(node.key), fontname="helvetica",shape="circle",color="red")
    for k in range(len(node.successor)):
        if node.successor[k]:
            label1=str([node.successor[k][0],node.successor[k][1]])
            if not node.successor[k][2] in succ:
                dot=layout(node.successor[k][2],dot,succ)
                dot.edge(str(node.idx),str(node.successor[k][2].idx),color=col[k%4],label=label1)
                succ.append(node.successor[k][2])
            else:
                dot.edge(str(node.idx),str(node.successor[k][2].idx),color=col[k%4],label=label1)
    return dot        

        
def Ini_BDD(index_order=[]):
    """To initialize the unique_table,computed_table and set up a global index order"""
    global computed_table
    global unique_table
    global global_node_idx
    global add_find_time,add_hit_time,cont_find_time,cont_hit_time
    global_node_idx=0

    unique_table = dict()
    computed_table = dict()
    add_find_time=0
    add_hit_time=0
    cont_find_time=0
    cont_hit_time=0
    set_index_order(index_order)

def Clear_BDD():
    """To initialize the unique_table,computed_table and set up a global index order"""
    global computed_table
    global unique_table
    global global_node_idx
    global add_find_time,add_hit_time,cont_find_time,cont_hit_time
    global_node_idx=0
    unique_table.clear()
    computed_table.clear()
    add_find_time=0
    add_hit_time=0
    cont_find_time=0
    cont_hit_time=0
    global_node_idx=0

def get_one_state():
    node = Find_Or_Add_Unique_table(-1)
    tdd = BDD(node)
    return tdd
def get_zero_state():
    node = Find_Or_Add_Unique_table(-1)
    tdd = BDD(node)
    tdd.weight=0
    return tdd

def get_unique_table():
    return unique_table

def get_unique_table_num():
    return len(unique_table)

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

def get_node_set(node,node_set=set()):
    """Only been used when counting the node number of a TDD"""
    if not node.idx in node_set:
        node_set.add(node.idx)
        for k in node.successor:
            node_set = get_node_set(k[2],node_set)
    return node_set

def Find_Or_Add_Unique_table(x,the_successors=[]):
    """To return a node if it already exist, creates a new node otherwise"""
    global global_node_idx,unique_table
    
    if x==-1:
        if unique_table.__contains__(x):
            return unique_table[x]
        else:
            res=Node(x)
            res.idx=0
            unique_table[x]=res
        return res
    temp_key=[x]
    for succ in the_successors:
        temp_key.append(succ[0])
        temp_key.append(get_int_key(succ[1]))
        temp_key.append(succ[2])
    temp_key=tuple(temp_key)
    if temp_key in unique_table:
        return unique_table[temp_key]
    else:
        res=Node(x)
        global_node_idx+=1
        res.idx=global_node_idx
        res.successor=the_successors
        
        res.index_degree={x:the_successors[0][0]}
        
        common_keys=the_successors[0][2].index_degree.keys()
        union_keys=the_successors[0][2].index_degree.keys()
        for succ in the_successors[1:]:
            common_keys=common_keys&succ[2].index_degree.keys()
            union_keys=union_keys|succ[2].index_degree.keys()
        
        for k in union_keys:
            res.index_degree[k]=0
        
        for k in common_keys:
            min_degree=the_successors[0][2].index_degree[k]
            for succ in the_successors[1:]:
                min_degree=min(succ[2].index_degree[k],min_degree)
            res.index_degree[k]=min_degree
        unique_table[temp_key]=res
    return res


def normalize(x,the_successors):
    """The normalize and reduce procedure"""
    global epi

    all_zero=True
    remove_list=[]
    for succ in the_successors:
        if get_int_key(succ[1])!=(0,0):
            all_zero=False
            # break
        else:
            remove_list.append(succ)

    for succ in remove_list:
        the_successors.remove(succ)

    if all_zero:
        res=BDD(Find_Or_Add_Unique_table(-1))
        res.weight=0
        return res

    if len(the_successors)==1 and the_successors[0][0]==0:
        if not get_int_key(the_successors[0][1])==(0,0):
            res=BDD(the_successors[0][2])
            res.weight=the_successors[0][1]
            # res.weight=nsimplify(the_successors[0][1],tolerance=1e-3)
        else:
            res=BDD(Find_Or_Add_Unique_table(-1))
            res.weight=0
        return res
    
    the_successors.sort()
    
    the_degrees=[succ[0] for succ in the_successors]
    if not len(set(the_degrees))==len(the_degrees):
        print('Repeated degrees')
    
    weigs_abs=[abs(succ[1])/epi for succ in the_successors]
    weig_max=the_successors[weigs_abs.index(max(weigs_abs))][1]
    
    the_successors=[[succ[0],succ[1]/weig_max,succ[2]] for succ in the_successors]
    
    node=Find_Or_Add_Unique_table(x,the_successors)
    res=BDD(node)
    res.weight=weig_max
    
    return res

def get_count():
    global add_find_time,add_hit_time,cont_find_time,cont_hit_time
    print("add:",add_hit_time,'/',add_find_time,'/',add_hit_time/add_find_time)
    print("cont:",cont_hit_time,"/",cont_find_time,"/",cont_hit_time/cont_find_time)

def find_computed_table(item):
    """To return the results that already exist"""
    global computed_table,add_find_time,add_hit_time,cont_find_time,cont_hit_time
    if item[0]=='s':
        temp_key=item[1].index_2_key[item[2]]
        the_key=('s',get_int_key(item[1].weight),item[1].node,temp_key,item[3])
        if computed_table.__contains__(the_key):
            res = computed_table[the_key]
            tdd = BDD(res[1])
            tdd.weight = res[0]
            return tdd
    elif item[0] == '+':
        the_key=('+',get_int_key(item[1].weight),item[1].node,get_int_key(item[2].weight),item[2].node)
        add_find_time+=1
        if computed_table.__contains__(the_key):
            res = computed_table[the_key]
            tdd = BDD(res[1])
            tdd.weight = res[0]
            add_hit_time+=1
            return tdd
        the_key=('+',get_int_key(item[2].weight),item[2].node,get_int_key(item[1].weight),item[1].node)
        if computed_table.__contains__(the_key):
            res = computed_table[the_key]
            tdd = BDD(res[1])
            tdd.weight = res[0]
            add_hit_time+=1
            return tdd
    elif item[0] == '*':
        the_key=('*',get_int_key(item[1].weight),item[1].node,get_int_key(item[2].weight),item[2].node)
        # ,*item[3]
        cont_find_time+=1
        if computed_table.__contains__(the_key):
            res = computed_table[the_key]
            tdd = BDD(res[1])
            tdd.weight = res[0]
            cont_hit_time+=1            
            return tdd
        the_key=('*',get_int_key(item[2].weight),item[2].node,get_int_key(item[1].weight),item[1].node)
        if computed_table.__contains__(the_key):
            res = computed_table[the_key]
            tdd = BDD(res[1])
            tdd.weight = res[0]
            cont_hit_time+=1            
            return tdd
    elif item[0] == '/':
        the_key = ('/',get_int_key(item[1].weight),item[1].node,get_int_key(item[2].weight),item[2].node)
        if computed_table.__contains__(the_key):
            res = computed_table[the_key]
            a = BDD(res[1])
            a.weight = res[0]
            b = BDD(res[3])
            b.weight = res[2]
            c = BDD(res[5])
            c.weight = res[4]
            return [a,b,c]
    elif item[0] == 'r':
        the_key = ('r',get_int_key(item[1].weight),item[1].node,(tuple(item[2].keys()),tuple(item[2].values())))
        if computed_table.__contains__(the_key):
            res = computed_table[the_key]
            tdd = BDD(res[1])
            tdd.weight = res[0]          
            return tdd            
    return None

def insert_2_computed_table(item,res):
    """To insert an item to the computed table"""
    global computed_table,cont_time,find_time,hit_time
    if item[0]=='s':
        temp_key=item[1].index_2_key[item[2]]
        the_key = ('s',get_int_key(item[1].weight),item[1].node,temp_key,item[3])
        computed_table[the_key] = (res.weight,res.node)
    elif item[0] == '+':
        the_key = ('+',get_int_key(item[1].weight),item[1].node,get_int_key(item[2].weight),item[2].node)
        computed_table[the_key] = (res.weight,res.node)
    elif item[0] == '*':
        the_key = ('*',get_int_key(item[1].weight),item[1].node,get_int_key(item[2].weight),item[2].node)
        computed_table[the_key] = (res.weight,res.node)
    elif item[0] == '/':
        the_key = ('/',get_int_key(item[1].weight),item[1].node,get_int_key(item[2].weight),item[2].node)
        computed_table[the_key] = (res[0].weight,res[0].node,res[1].weight,res[1].node,res[2].weight,res[2].node)
    elif item[0] == 'r':
        the_key = ('r',get_int_key(item[1].weight),item[1].node,(tuple(item[2].keys()),tuple(item[2].values())))
        computed_table[the_key] = (res.weight,res.node)

def get_bdd(function):
    global global_index_order 
    if isinstance(function,int) or isinstance(function,float) or isinstance(function,complex):
        bdd=get_one_state()
        bdd.weight=function
        return bdd
    try:
        function.args
    except:
        bdd=get_one_state()
        bdd.weight=complex(function)
        return bdd
    if len(function.args)==0:
        bdd=get_one_state()
        bdd.weight=complex(function)
        return bdd
#不能使用str(function)
    if len(function.args)==1:
        '''目前假設只有一個變數，如果變多要改'''
        
        bdd=normalize(str(function),[[1,1,Find_Or_Add_Unique_table(-1)]])
        
        #在這把symbol次序放入
        for item in function.free_symbols:
            if '[' in str(function.free_symbols):
                if '[' in str(item):
                    symbol_name=item.name
            else:
                symbol_name=item.name

        bdd.key_2_index['sin(%s)'%symbol_name]=global_index_order['sin(%s)'%symbol_name]
        bdd.key_2_index['cos(%s)'%symbol_name]=global_index_order['cos(%s)'%symbol_name]
        print('BDD 463',bdd.key_2_index)
        return bdd
    if isinstance(function,sp.core.add.Add):
        bdd=get_zero_state()
        for add_term in function.args:
            temp_tdd=get_bdd(add_term)
            bdd=add(bdd,temp_tdd)
    elif isinstance(function,sp.core.mul.Mul):
    # elif isinstance(function,sp.mul.Mul):
        bdd=get_one_state()
        for mul_term in function.args:
            temp_tdd=get_bdd(mul_term)
            bdd=mul(bdd,temp_tdd)
    elif isinstance(function,sp.core.power.Pow):
    # elif isinstance(function,sp.power.Pow):
        bdd=get_one_state()
        pow= function.args[1]
        for mul_times in range(pow):
            temp_tdd=get_bdd(function.args[0])
            bdd=mul(bdd,temp_tdd)
    return bdd



def mul(tdd1,tdd2):
    """The contraction of two TDDs, var_cont is in the form [[4,1],[3,2]]"""
    global global_index_order 

    k1=tdd1.node.key
    k2=tdd2.node.key
    w1=tdd1.weight
    w2=tdd2.weight
    
    if k1==k2==-1:
        weig=tdd1.weight*tdd2.weight
        term=Find_Or_Add_Unique_table(-1)
        res=BDD(term)        
        if get_int_key(weig)==(0,0):
            res.weight=0
            return res
        else:
            res.weight=weig
            return res        

    if k1==-1:
        if w1==0:
            tdd=BDD(tdd1.node)
            tdd.weight=0
            return tdd

        tdd=BDD(tdd2.node)
        tdd.weight=w1*w2
        return tdd
            
    if k2==-1:
        if w2==0:
            tdd=BDD(tdd2.node)
            tdd.weight=0
            return tdd        

        tdd=BDD(tdd1.node)
        tdd.weight=w1*w2
        return tdd
    
    tdd1.weight=1
    tdd2.weight=1
    
    tdd=find_computed_table(['*',tdd1,tdd2])
    if tdd:
        tdd.weight=tdd.weight*w1*w2
        tdd1.weight=w1
        tdd2.weight=w2        
        return tdd
              
    tdd=BDD(tdd2.node)
    tdd.weight=0
    if k1==k2:
        for succ1 in tdd1.node.successor:
            for succ2 in tdd2.node.successor:
                temp_tdd1=BDD(succ1[2])
                temp_tdd1.weight=succ1[1]
                temp_tdd2=BDD(succ2[2])
                temp_tdd2.weight=succ2[1]            
                temp_res=mul(temp_tdd1,temp_tdd2)
                if not succ1[0]+succ2[0]==0:
                    if succ1[0]+succ2[0]==2 and k1[:3]=='sin':
                        temp_res1=normalize('cos'+k1[3:],[[2,-1,Find_Or_Add_Unique_table(-1)]])
                        temp_res1=mul(temp_res1,temp_res)
                        temp_res=add(temp_res,temp_res1)
                    else:
                        temp_res=normalize(k1,[[succ1[0]+succ2[0],temp_res.weight,temp_res.node]])
                tdd=add(tdd,temp_res)
    elif global_index_order[k1]<=global_index_order[k2]:
        the_successor=[]
        for succ1 in tdd1.node.successor:
            temp_tdd1=BDD(succ1[2])
            temp_tdd1.weight=succ1[1]
            temp_res=mul(temp_tdd1,tdd2)
            the_successor.append([succ1[0],temp_res.weight,temp_res.node])
            tdd=normalize(k1,the_successor)
    else:
        the_successor=[]
        for succ2 in tdd2.node.successor:
            temp_tdd2=BDD(succ2[2])
            temp_tdd2.weight=succ2[1]
            temp_res=mul(tdd1,temp_tdd2)
            the_successor.append([succ2[0],temp_res.weight,temp_res.node])
            tdd=normalize(k2,the_successor)        
            
    insert_2_computed_table(['*',tdd1,tdd2],tdd)
    tdd.weight=tdd.weight*w1*w2
    tdd1.weight=w1
    tdd2.weight=w2
    
    return tdd      
        

def add(tdd1,tdd2):
    """The apply function of two TDDs. Mostly, it is used to do addition here."""
    global global_index_order    

    if tdd1.weight==0:
        return tdd2.self_copy()
    
    if tdd2.weight==0:
        return tdd1.self_copy()
    
    if tdd1.node==tdd2.node:
        weig=tdd1.weight+tdd2.weight
        if get_int_key(weig)==(0,0):
            term=Find_Or_Add_Unique_table(-1)
            res=BDD(term)
            res.weight=0
            return res
        else:
            res=BDD(tdd1.node)
            res.weight=weig
            return res

    if find_computed_table(['+',tdd1,tdd2]):
        return find_computed_table(['+',tdd1,tdd2])
    
    k1=tdd1.node.key
    k2=tdd2.node.key
    
    the_successor=[]
    if k1==k2:
        the_key=k1
        degree1=[succ[0] for succ in tdd1.node.successor]
        degree2=[succ[0] for succ in tdd2.node.successor]
        the_successor+=[[succ[0],succ[1]*tdd1.weight,succ[2]] for succ in tdd1.node.successor if not succ[0] in degree2]
        the_successor+=[[succ[0],succ[1]*tdd2.weight,succ[2]] for succ in tdd2.node.successor if not succ[0] in degree1]
        com_degree=[d for d in degree1 if d in degree2]
        for d in com_degree:
            pos1=degree1.index(d)
            temp_tdd1=BDD(tdd1.node.successor[pos1][2])
            temp_tdd1.weight=tdd1.node.successor[pos1][1]*tdd1.weight
            pos2=degree2.index(d)
            temp_tdd2=BDD(tdd2.node.successor[pos2][2])
            temp_tdd2.weight=tdd2.node.successor[pos2][1]*tdd2.weight          
            temp_res=add(temp_tdd1,temp_tdd2)
#             print(temp_tdd1,temp_tdd2,temp_res)
            the_successor.append([d,temp_res.weight,temp_res.node])
    elif global_index_order[k1]<=global_index_order[k2]:
        the_key=k1
        if tdd1.node.successor[0][0]==0:
            temp_tdd1=BDD(tdd1.node.successor[0][2])
            temp_tdd1.weight=tdd1.node.successor[0][1]*tdd1.weight
            temp_res=add(temp_tdd1,tdd2)
            the_successor.append([0,temp_res.weight,temp_res.node])
            the_successor+=[[succ[0],succ[1]*tdd1.weight,succ[2]] for succ in tdd1.node.successor[1:]]
        else:
            the_successor.append([0,tdd2.weight,tdd2.node])
            the_successor+=[[succ[0],succ[1]*tdd1.weight,succ[2]] for succ in tdd1.node.successor]
    else:
        the_key=k2
        if tdd2.node.successor[0][0]==0:
            temp_tdd2=BDD(tdd2.node.successor[0][2])
            temp_tdd2.weight=tdd2.node.successor[0][1]*tdd2.weight
            temp_res=add(tdd1,temp_tdd2)
            the_successor.append([0,temp_res.weight,temp_res.node])
            the_successor+=[[succ[0],succ[1]*tdd2.weight,succ[2]] for succ in tdd2.node.successor[1:]]
        else:
            the_successor.append([0,tdd1.weight,tdd1.node])
            the_successor+=[[succ[0],succ[1]*tdd2.weight,succ[2]] for succ in tdd2.node.successor]          
            
    res = normalize(the_key,the_successor)
    insert_2_computed_table(['+',tdd1,tdd2],res)
#     print(res)
#     print('-----------')
    return res

def reduce_degree(tdd,min_degree):
    if len(min_degree)==0:
        return tdd.self_copy()
    
    max_index=max([global_index_order[k] for k in min_degree])
    
    k1=tdd.node.key
    
    if global_index_order[k1]>max_index:
        return tdd.self_copy()
    
    w=tdd.weight
    tdd.weight=1
    
    if find_computed_table(['r',tdd,min_degree]):
        res=find_computed_table(['r',tdd,min_degree])
        res.weight*=w
        tdd.weight=w
        return res
    
    the_successor=[]
    
    if k1 in min_degree:
        for succ in tdd.node.successor:
            temp_res=reduce_degree(BDD(succ[2]),min_degree)
            the_successor.append([succ[0]-min_degree[k1],succ[1]*temp_res.weight,temp_res.node])
        res=normalize(k1,the_successor)
    else:
        for succ in tdd.node.successor:
            temp_res=reduce_degree(BDD(succ[2]),min_degree)
            the_successor.append([succ[0],succ[1]*temp_res.weight,temp_res.node])
        res=normalize(k1,the_successor)
    insert_2_computed_table(['r',tdd,min_degree],res)
    res.weight*=w
    tdd.weight=w    
    return res


def normalize_2_fun(tdd1,tdd2):
    #tdd2/tdd1,the result is like [tdd1,1,tdd2/tdd1]
    if tdd1.weight==0:
        return [tdd2.self_copy(),tdd1,get_one_state()]
    
    if tdd2.weight==0:
        return [tdd1.self_copy(),get_one_state(),tdd2.self_copy()]    
    
    if tdd1.node.key==-1:
        temp=tdd2.self_copy()
        temp.weight/=tdd1.weight
        return [tdd1.self_copy(),get_one_state(),temp]
    
    if tdd1.node==tdd2.node:
        temp=get_one_state()
        temp.weight=tdd2.weight/tdd1.weight
        return [tdd1.self_copy(),get_one_state(),temp]
    
    # return [get_one_state(),tdd1,tdd2]

    w1=tdd1.weight
    w2=tdd2.weight
    tdd1.weight=1
    tdd2.weight=1
    
    if find_computed_table(['/',tdd1,tdd2]):
        [a,b,c]=find_computed_table(['/',tdd1,tdd2])
        a.weight*=w1
        c.weight*=w2/w1
        tdd1.weight=w1
        tdd2.weight=w2
        return [a,b,c]
    
    min_degree=dict()
    for key in tdd1.node.index_degree:
        if key in tdd2.node.index_degree:
            if min(tdd1.node.index_degree[key],tdd2.node.index_degree[key])>0:
                min_degree[key]=min(tdd1.node.index_degree[key],tdd2.node.index_degree[key])
                
    a=get_one_state()                
    if len(min_degree)>0:
#         print(min_degree)
        b=reduce_degree(tdd1,min_degree)
        c=reduce_degree(tdd2,min_degree)
        min_degree_order=[[global_index_order[k],k] for k in min_degree]
        while len(min_degree_order)>0:
            key=min_degree_order.pop()[1]
            a=normalize(key,[[min_degree[key],1,a.node]])        
    else:
        b=tdd1.self_copy()
        c=tdd2.self_copy()
    
    a.weight=b.weight
    b.weight=1
    c.weight/=a.weight
#     w=max(tdd1.weight,tdd2.weight)
#     if not get_int_key(tdd1.weight)==(0,0):
#         w=tdd1.weight
#     else:
#         w=tdd2.weight
#     a=BDD(Find_Or_Add_Unique_table(-1))
#     a.weight=w
#     b=tdd1.self_copy()
#     b.weight=tdd1.weight/w
#     c=tdd2.self_copy()
#     c.weight=tdd2.weight/w
    insert_2_computed_table(['/',tdd1,tdd2],[a,b,c])
    
    a.weight*=w1
    c.weight*=w2/w1
    tdd1.weight=w1
    tdd2.weight=w2
    return [a,b,c]        
    
    


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
        # print(734,expr)
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


def var_sort (var):
    global global_index_order, inverse_global_index_order
    
    idx=[global_index_order[k] for k in var]
    idx.sort( )
    var_sort= [inverse_global_index_order[k] for k in idx]

    return var_sort

def cont(bdd1,bdd2):
    #找出哪些要輸出(cont)/保留(out)
    var_out1=[var for var in bdd1.index_set ]
    var_out2=[var for var in bdd2.index_set ]

    var_out=var_out1+var_out2
    var_out=var_sort(var_out) #Index 已含有比較大小的函數

    # var_out_idx=[var.key for var in var_out]
    # var_cont_idx=[var.key for var in var_cont]
    # var_cont_idx=[var for var in var_cont_idx if not var in var_out_idx]
    
    idx_2_key={-1:-1}
    key_2_idx={-1:-1}
    
    n=0
    for k in range(len(var_out)-1,-1,-1):
        if not var_out[k] in idx_2_key:
            idx_2_key[var_out[k]]=n
            key_2_idx[n]=var_out[k]
            n+=1
    #找key的對應
    key_2_new_key=[[],[]]
    #找cont的順序
    cont_order=[[],[]]
    for k in range(len(bdd1.key_2_index)-1):
        v=bdd1.key_2_index[k]
        if v in idx_2_key:
            key_2_new_key[0].append(idx_2_key[v])
        else:
            key_2_new_key[0].append('c')
        cont_order[0].append(global_index_order[v])
        
    cont_order[0].append(float('inf'))
    
    for k in range(len(bdd2.key_2_index)-1):     
        v=bdd2.key_2_index[k]
        if v in idx_2_key:
            key_2_new_key[1].append(idx_2_key[v])
        else:
            key_2_new_key[1].append('c')
        cont_order[1].append(global_index_order[v])
    cont_order[1].append(float('inf'))

    tdd=mul2(bdd1,bdd2,key_2_new_key,cont_order)
    tdd.index_set=var_out
    tdd.index_2_key=idx_2_key
    tdd.key_2_index=key_2_idx

    return tdd


def mul2(bdd1,bdd2,key_2_new_key,cont_order):
    """The contraction of two TDDs, var_cont is in the form [[4,1],[3,2]]"""
    global global_index_order ,inverse_global_index_order

    k1=bdd1.node.key
    k2=bdd2.node.key
    w1=bdd1.weight
    w2=bdd2.weight
    
    if k1==k2==-1:
        weig=bdd1.weight*bdd2.weight
        term=Find_Or_Add_Unique_table(-1)
        res=BDD(term)        
        if get_int_key(weig)==(0,0):
            res.weight=0
            return res
        else:
            res.weight=weig
            return res        

    if k1==-1:
        if w1==0:
            tdd=BDD(bdd1.node)
            tdd.weight=0
            return tdd

        tdd=BDD(bdd2.node)
        tdd.weight=w1*w2
        return tdd
            
    if k2==-1:
        if w2==0:
            tdd=BDD(bdd2.node)
            tdd.weight=0
            return tdd        

        tdd=BDD(bdd1.node)
        tdd.weight=w1*w2
        return tdd
    
    bdd1.weight=1
    bdd2.weight=1
    
    temp_key_2_new_key=[]
    temp_key_2_new_key.append(tuple([k for k in key_2_new_key[0][:(k1+1)]]))
    temp_key_2_new_key.append(tuple([k for k in key_2_new_key[1][:(k2+1)]]))

    tdd=find_computed_table(['*',bdd1,bdd2,temp_key_2_new_key])
    if tdd:
        tdd.weight=tdd.weight*w1*w2
        bdd1.weight=w1
        bdd2.weight=w2        
        return tdd
              
    tdd=BDD(bdd2.node)
    tdd.weight=0
    if cont_order[0][k1] == cont_order[1][k2]:
        for succ1 in bdd1.node.successor:
            for succ2 in bdd2.node.successor:
                temp_bdd1=BDD(succ1[2])
                temp_bdd1.weight=succ1[1]
                temp_tdd2=BDD(succ2[2])
                temp_tdd2.weight=succ2[1]            
                temp_res=mul2(temp_bdd1,temp_tdd2,key_2_new_key,cont_order)
                if not succ1[0]+succ2[0]==0:
                    if succ1[0]+succ2[0]==2 and inverse_global_index_order[cont_order[0][k1]][:3]=='sin': # 這邊還需要改
                        temp_res1=normalize(key_2_new_key[0][k1]-1,[[2,-1,Find_Or_Add_Unique_table(-1)]])
                        temp_res1=mul2(temp_res1,temp_res,key_2_new_key,cont_order)
                        temp_res=add2(temp_res,temp_res1,key_2_new_key,cont_order)
                    else:
                        temp_res=normalize(key_2_new_key[0][k1],[[succ1[0]+succ2[0],temp_res.weight,temp_res.node]])
                tdd=add2(tdd,temp_res,key_2_new_key,cont_order)
    elif cont_order[0][k1] <= cont_order[1][k2]:
        the_successor=[]
        for succ1 in bdd1.node.successor:
            temp_bdd1=BDD(succ1[2])
            temp_bdd1.weight=succ1[1]
            temp_res=mul2(temp_bdd1, bdd2, key_2_new_key, cont_order)
            the_successor.append([succ1[0],temp_res.weight,temp_res.node])
            tdd=normalize(key_2_new_key[0][k1],the_successor)
    else:
        the_successor=[]
        for succ2 in bdd2.node.successor:
            temp_tdd2=BDD(succ2[2])
            temp_tdd2.weight=succ2[1]
            temp_res=mul2(bdd1, temp_tdd2 ,key_2_new_key, cont_order)
            the_successor.append([succ2[0],temp_res.weight,temp_res.node])
            tdd=normalize(key_2_new_key[1][k2],the_successor)        
            
    insert_2_computed_table(['*',bdd1,bdd2],tdd)
    tdd.weight=tdd.weight*w1*w2
    bdd1.weight=w1
    bdd2.weight=w2
    
    return tdd      


def add2(bdd1,bdd2,key_2_new_key,cont_order):
    """The apply function of two TDDs. Mostly, it is used to do addition here."""
    global global_index_order    

    if bdd1.weight==0:
        return bdd2.self_copy()
    
    if bdd2.weight==0:
        return bdd1.self_copy()
    
    if bdd1.node==bdd2.node:
        weig=bdd1.weight+bdd2.weight
        if get_int_key(weig)==(0,0):
            term=Find_Or_Add_Unique_table(-1)
            res=BDD(term)
            res.weight=0
            return res
        else:
            res=BDD(bdd1.node)
            res.weight=weig
            return res

    if find_computed_table(['+',bdd1,bdd2]):
        return find_computed_table(['+',bdd1,bdd2])
    
    k1=bdd1.node.key
    k2=bdd2.node.key
    
    the_successor=[]
    if k1==k2:
        the_key=k1
        degree1=[succ[0] for succ in bdd1.node.successor]
        degree2=[succ[0] for succ in bdd2.node.successor]
        the_successor+=[[succ[0],succ[1]*bdd1.weight,succ[2]] for succ in bdd1.node.successor if not succ[0] in degree2]
        the_successor+=[[succ[0],succ[1]*bdd2.weight,succ[2]] for succ in bdd2.node.successor if not succ[0] in degree1]
        com_degree=[d for d in degree1 if d in degree2]
        for d in com_degree:
            pos1=degree1.index(d)
            temp_tdd1=BDD(bdd1.node.successor[pos1][2])
            temp_tdd1.weight=bdd1.node.successor[pos1][1]*bdd1.weight
            pos2=degree2.index(d)
            temp_tdd2=BDD(bdd2.node.successor[pos2][2])
            temp_tdd2.weight=bdd2.node.successor[pos2][1]*bdd2.weight          
            temp_res=add2(temp_tdd1,temp_tdd2,key_2_new_key,cont_order)
#             print(temp_tdd1,temp_tdd2,temp_res)
            the_successor.append([d,temp_res.weight,temp_res.node])
    elif global_index_order[k1]<=global_index_order[k2]:
        the_key=k1
        if bdd1.node.successor[0][0]==0:
            temp_tdd1=BDD(bdd1.node.successor[0][2])
            temp_tdd1.weight=bdd1.node.successor[0][1]*bdd1.weight
            temp_res=add2(temp_tdd1,bdd2,key_2_new_key,cont_order)
            the_successor.append([0,temp_res.weight,temp_res.node])
            the_successor+=[[succ[0],succ[1]*bdd1.weight,succ[2]] for succ in bdd1.node.successor[1:]]
        else:
            the_successor.append([0,bdd2.weight,bdd2.node])
            the_successor+=[[succ[0],succ[1]*bdd1.weight,succ[2]] for succ in bdd1.node.successor]
    else:
        the_key=k2
        if bdd2.node.successor[0][0]==0:
            temp_tdd2=BDD(bdd2.node.successor[0][2])
            temp_tdd2.weight=bdd2.node.successor[0][1]*bdd2.weight
            temp_res=add2(bdd1,temp_tdd2,key_2_new_key,cont_order)
            the_successor.append([0,temp_res.weight,temp_res.node])
            the_successor+=[[succ[0],succ[1]*bdd2.weight,succ[2]] for succ in bdd2.node.successor[1:]]
        else:
            the_successor.append([0,bdd1.weight,bdd1.node])
            the_successor+=[[succ[0],succ[1]*bdd2.weight,succ[2]] for succ in bdd2.node.successor]          
            
    res = normalize(the_key,the_successor)
    insert_2_computed_table(['+',bdd1,bdd2],res)
#     print(res)
#     print('-----------')
    return res