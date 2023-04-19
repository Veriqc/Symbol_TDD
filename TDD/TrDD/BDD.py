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
inverse_global_index_order = dict()
global_node_idx=0
add_find_time=0
add_hit_time=0
cont_find_time=0
cont_hit_time=0
epi=1e-10

    
class Node:
    """To define the node of bdd"""
    def __init__(self,key):
        self.idx = 0
        self._key = key
        self.successor = [] #the element should be in this form [degree, weight, node]
        self.index_degree = dict()
        self.value = dict()
        self.succ_num = 2
        self._repr = None

    def __repr__(self):
        if self._repr:
            return self._repr
        self._repr = str(self.key)+str(self.successor)
        return self._repr
    def __hash__(self) -> int:
        return hash(self.__repr__())
    
    @property 
    def key(self):
        if not isinstance(self._key, int):
            print('BDD 43', self._key)
            self._key = int(self._key.real)
        return self._key

    def expr(self, key_2_index):
        global inverse_global_index_order

        if self.key==-1:
            return 1+0j

        def get_numbers(s):
            import re
            # 使用正則表達式找到所有匹配"\d+"的子串，即連續的一個或多個數字
            numbers = re.findall("\d+", s)
            return int(numbers[0])
        # param_expr=inverse_global_index_order[self.key]
        param_expr = key_2_index[self.key]
        sym_str=param_expr.replace('sin(','').replace('cos(','').replace('[','').replace('])','').replace(str(get_numbers(param_expr)),'')

        sym_base = IndexedBase(sym_str)

        sp_expr=sympify(str(param_expr), locals={sym_str: sym_base})
        
        res=0
        for succ in self.successor:
            res+=nsimplify(succ[1]*sp_expr**succ[0],tolerance=1e-3)*succ[2].expr(key_2_index)
        self._expr=res
        return res
    
    


class BDD:
    def __init__(self,node, weight=1+0j ,key_2_index= {-1:-1}):
        """BDD"""
        self._weight = weight
        self._expr=None
        self.key_2_index = key_2_index

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
        temp.key_2_index = self.key_2_index.copy()
        temp._expr = self._expr
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
        print('BDD 129', str(self.expr) , str(other.expr))
        return str(self.expr) == str(other.expr)
    # self.node == other.node and \
    #         math.isclose(self.weight.real, other.weight.real, rel_tol=epi) and \
    #         math.isclose(self.weight.imag, other.weight.imag, rel_tol=epi) and \
    #         self.key_2_index == other.key_2_index

    def __add__(self, g):        
        # return add(self,g)
         return cont('add',self, g)
        # return add3(self,g)
    
    def __mul__(self, g):
        # return mul(self,g)
        return cont('mul',self, g)
        # return mul3(self, g)
    
    def self_normalize(self, g):
        return normalize_2_fun(g,self)
    
    @property
    def expr(self):
        if self._expr:
            return self._expr

        value= np.round(self.weight,int(-np.log10(epi)))

        if self.node.key==-1:
            return value
        
        res=self.node.expr(self.key_2_index)
        self._expr = value*res
        return self._expr
    
    def get_value(self,val):
        res=get_value_node(self.node,val)
        return self.weight*res
    
    def get_value2(self,val):
        res=get_value_node2(self.node,val)
        return self.weight*res    
        
    def get_key_2_index_subdict(self, top_node, subdict=dict()): 
        if top_node.key >=0 :
            cos_key = (top_node.key//2)*2
            sin_key =  cos_key+1
            subdict[cos_key] = self.key_2_index[cos_key]
            subdict[sin_key] = self.key_2_index[sin_key]
        else:
            subdict[top_node.key] = self.key_2_index[top_node.key]
        for succ in top_node.successor:
            self.get_key_2_index_subdict(succ[2], subdict)
        return subdict

    def __hash__(self):
        return hash(self.__repr__())
    
    def __repr__(self):
        return str(self.expr)
        
    
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
    global global_index_order, inverse_global_index_order
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
    bdd = BDD(node)
    bdd.key_2_index = {-1:-1}
    return bdd
    # return BDD(1)
def get_zero_state():
    node = Find_Or_Add_Unique_table(-1)
    bdd = BDD(node)
    bdd.weight=0
    bdd.key_2_index = {-1:-1}
    return bdd
    # return BDD(0)
def get_unique_table():
    return unique_table

def get_unique_table_num():
    return len(unique_table)

def set_index_order(var_order):
    global global_index_order, inverse_global_index_order
    # index: sin, cos, ; order: the cont order 0,1 ..., inf
    global_index_order=dict()
    if isinstance(var_order,list):
        for k in range(len(var_order)):
            global_index_order[var_order[k]]=k
    if isinstance(var_order,dict):
        global_index_order = copy.copy(var_order)
    global_index_order[-1] = float('inf')

    inverse_global_index_order = {v: k for k, v in global_index_order.items()}
    
def get_index_order():
    global global_index_order
    return copy.copy(global_index_order)
    
def get_int_key(weight):
    """To transform a complex number to a tuple with int values"""
    global epi

    return (int(round(weight.real/epi)) ,int(round(weight.imag/epi)))

def get_node_set(node,node_set=set()):
    """Only been used when counting the node number of a bdd"""
    if not node.idx in node_set:
        node_set.add(node.idx)
        for k in node.successor:
            node_set = get_node_set(k[2],node_set)
    return node_set

def Find_Or_Add_Unique_table(x,the_successors=[]):
    """To return a node if it already exist, creates a new node otherwise"""
    global global_node_idx,unique_table
    if not isinstance(x,int):
        print('BDD 284',x)
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


def normalize(key, the_successors):
    """The normalize and reduce procedure"""
    if not isinstance(key, int):
        print('BDD 330', key)
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
        return BDD(Find_Or_Add_Unique_table(-1),weight=0)

    if len(the_successors)==1 and the_successors[0][0]==0:
        if not get_int_key(the_successors[0][1])==(0,0):
            succ = the_successors[0]
            res=BDD(succ[2], weight=succ[1])
        else:
            res = BDD(Find_Or_Add_Unique_table(-1),weight=0)
        return res
    
    the_successors.sort()
    
    the_degrees=[succ[0] for succ in the_successors]
    if not len(set(the_degrees))==len(the_degrees):
        print('Repeated degrees')
    
    weigs_abs=[abs(succ[1])/epi for succ in the_successors]
    weig_max=the_successors[weigs_abs.index(max(weigs_abs))][1]
    
    the_successors=[[succ[0],succ[1]/weig_max,succ[2]] for succ in the_successors]
    
    node=Find_Or_Add_Unique_table(key, the_successors)
    res=BDD(node ,weight=weig_max)
    
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
            bdd = BDD(res[1])
            bdd.weight = res[0]
            return bdd
    elif item[0] == '+':
        the_key=('+',get_int_key(item[1].weight),item[1].node,get_int_key(item[2].weight),item[2].node,item[3][0],item[3][1])
        add_find_time+=1
        if computed_table.__contains__(the_key):
            res = computed_table[the_key]
            bdd = BDD(res[1])
            bdd.weight = res[0]
            add_hit_time+=1
            return bdd
        the_key=('+',get_int_key(item[2].weight),item[2].node,get_int_key(item[1].weight),item[1].node,item[3][0],item[3][1])
        if computed_table.__contains__(the_key):
            res = computed_table[the_key]
            bdd = BDD(res[1])
            bdd.weight = res[0]
            add_hit_time+=1
            return bdd
    elif item[0] == '*':
        the_key=('*',get_int_key(item[1].weight),item[1].node,get_int_key(item[2].weight),item[2].node,item[3][0],item[3][1])
        # ,*item[3]
        cont_find_time+=1
        if computed_table.__contains__(the_key):
            res = computed_table[the_key]
            bdd = BDD(res[1])
            bdd.weight = res[0]
            cont_hit_time+=1            
            return bdd
        the_key=('*',get_int_key(item[2].weight),item[2].node,get_int_key(item[1].weight),item[1].node,item[3][0],item[3][1])
        if computed_table.__contains__(the_key):
            res = computed_table[the_key]
            bdd = BDD(res[1])
            bdd.weight = res[0]
            cont_hit_time+=1            
            return bdd
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
            bdd = BDD(res[1])
            bdd.weight = res[0]          
            return bdd            
    return None

def insert_2_computed_table(item,res):
    """To insert an item to the computed table"""
    global computed_table,cont_time,find_time,hit_time
    if item[0]=='s':
        temp_key=item[1].index_2_key[item[2]]
        the_key = ('s',get_int_key(item[1].weight),item[1].node,temp_key,item[3])
        computed_table[the_key] = (res.weight,res.node)
    elif item[0] == '+':
        the_key = ('+',get_int_key(item[1].weight),item[1].node,get_int_key(item[2].weight),item[2].node,item[3][0],item[3][1])
        computed_table[the_key] = (res.weight,res.node)
    elif item[0] == '*':
        the_key = ('*',get_int_key(item[1].weight),item[1].node,get_int_key(item[2].weight),item[2].node,item[3][0],item[3][1])
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
        # print('BDD 465',function)
        
        #在這把symbol次序放入
        for item in function.free_symbols:
            if '[' in str(function.free_symbols):
                if '[' in str(item):
                    symbol_name=item.name
            else:
                symbol_name=item.name
        str_func = str(function)
        order = global_index_order[str_func]
        
        # bdd=normalize( order % 2, [[1,1,Find_Or_Add_Unique_table(-1)]])
        bdd = BDD(Node(order % 2))
        bdd.node.successor = [[1,1,Find_Or_Add_Unique_table(-1)]]
        bdd.key_2_index[0]='cos(%s)'%symbol_name
        bdd.key_2_index[1]='sin(%s)'%symbol_name
        bdd.key_2_index[-1]=-1

        print('BDD 507',bdd)
        return bdd
    # print('BDD 483', type(function))
    if isinstance(function,sp.core.add.Add):
        bdd=get_zero_state()
        for add_term in function.args:
            temp_bdd = get_bdd(add_term)
            bdd = bdd + temp_bdd
    elif isinstance(function,sp.core.mul.Mul):
    # elif isinstance(function,sp.mul.Mul):
        bdd=get_one_state()
        for mul_term in function.args:
            temp_bdd = get_bdd(mul_term)
            bdd = bdd * temp_bdd
    elif isinstance(function,sp.core.power.Pow):
    # elif isinstance(function,sp.power.Pow):
        bdd=get_one_state()
        pow= function.args[1]
        for mul_times in range(pow):
            temp_bdd = get_bdd(function.args[0])
            bdd = bdd * temp_bdd
    return bdd



def mul(bdd1,bdd2):
    """The contraction of two bdds, var_cont is in the form [[4,1],[3,2]]"""
    global global_index_order 

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
            bdd=BDD(bdd1.node)
            bdd.weight=0
            return bdd

        bdd=BDD(bdd2.node)
        bdd.weight=w1*w2
        return bdd
            
    if k2==-1:
        if w2==0:
            bdd=BDD(bdd2.node)
            bdd.weight=0
            return bdd        

        bdd=BDD(bdd1.node)
        bdd.weight=w1*w2
        return bdd
    
    bdd1.weight=1
    bdd2.weight=1
    
    bdd=find_computed_table(['*',bdd1,bdd2])
    if bdd:
        bdd.weight=bdd.weight*w1*w2
        bdd1.weight=w1
        bdd2.weight=w2        
        return bdd
              
    bdd=BDD(bdd2.node)
    bdd.weight=0
    if k1==k2:
        for succ1 in bdd1.node.successor:
            for succ2 in bdd2.node.successor:
                temp_bdd1=BDD(succ1[2])
                temp_bdd1.weight=succ1[1]
                temp_bdd2=BDD(succ2[2])
                temp_bdd2.weight=succ2[1]            
                temp_res=mul(temp_bdd1,temp_bdd2)
                if not succ1[0]+succ2[0]==0:
                    # if succ1[0]+succ2[0]==2 and k1[:3]=='sin':
                    if succ1[0]+succ2[0]==2 and k1 % 2 == 0 :
                        temp_res1=normalize( k1 + 1,[[2,-1,Find_Or_Add_Unique_table(-1)]])
                        temp_res1=mul(temp_res1,temp_res)
                        temp_res=add(temp_res,temp_res1)
                    else:
                        temp_res=normalize(k1,[[succ1[0]+succ2[0],temp_res.weight,temp_res.node]])
                bdd=add(bdd,temp_res)
    elif k1 <= k2:
    # elif global_index_order[k1]<=global_index_order[k2]:
        the_successor=[]
        for succ1 in bdd1.node.successor:
            temp_bdd1=BDD(succ1[2])
            temp_bdd1.weight=succ1[1]
            temp_res=mul(temp_bdd1,bdd2)
            the_successor.append([succ1[0],temp_res.weight,temp_res.node])
            bdd=normalize(k1,the_successor)
    else:
        the_successor=[]
        for succ2 in bdd2.node.successor:
            temp_bdd2=BDD(succ2[2])
            temp_bdd2.weight=succ2[1]
            temp_res=mul(bdd1,temp_bdd2)
            the_successor.append([succ2[0],temp_res.weight,temp_res.node])
            bdd=normalize(k2,the_successor)        
            
    insert_2_computed_table(['*',bdd1,bdd2],bdd)
    bdd.weight=bdd.weight*w1*w2
    bdd1.weight=w1
    bdd2.weight=w2
    
    return bdd      
        

def add(bdd1,bdd2):
    """The apply function of two bdds. Mostly, it is used to do addition here."""
    global global_index_order    

    if bdd1.weight == 0:
        return bdd2.self_copy()
    
    if bdd2.weight == 0:
        return bdd1.self_copy()
    
    if bdd1.node == bdd2.node:
        weig = bdd1.weight+bdd2.weight
        if get_int_key(weig) == (0,0):
            term = Find_Or_Add_Unique_table(-1)
            res = BDD(term)
            res.weight=0
            return res
        else:
            res = BDD(bdd1.node)
            res.weight = weig
            return res

    if find_computed_table(['+',bdd1,bdd2]):
        return find_computed_table(['+',bdd1,bdd2])
    
    k1 = bdd1.node.key
    k2 = bdd2.node.key
    print('BDD 627', k1, k2)
    the_successor = []
    if k1 == k2:
        the_key = k1
        degree1=[succ[0] for succ in bdd1.node.successor]
        degree2=[succ[0] for succ in bdd2.node.successor]
        the_successor+=[[succ[0],succ[1]*bdd1.weight,succ[2]] for succ in bdd1.node.successor if not succ[0] in degree2]
        the_successor+=[[succ[0],succ[1]*bdd2.weight,succ[2]] for succ in bdd2.node.successor if not succ[0] in degree1]
        com_degree=[d for d in degree1 if d in degree2]
        for d in com_degree:
            pos1=degree1.index(d)
            temp_bdd1=BDD(bdd1.node.successor[pos1][2])
            temp_bdd1.weight=bdd1.node.successor[pos1][1]*bdd1.weight
            pos2=degree2.index(d)
            temp_bdd2=BDD(bdd2.node.successor[pos2][2])
            temp_bdd2.weight=bdd2.node.successor[pos2][1]*bdd2.weight          
            temp_res=add(temp_bdd1,temp_bdd2)
#             print(temp_bdd1,temp_bdd2,temp_res)
            the_successor.append([d,temp_res.weight,temp_res.node])
    elif k1 <= k2:
    # elif global_index_order[k1]<=global_index_order[k2]:
        the_key = k1
        print('BDD 649',bdd1.node)
        if bdd1.node.successor[0][0] == 0:
            temp_bdd1=BDD(bdd1.node.successor[0][2])
            temp_bdd1.weight=bdd1.node.successor[0][1] * bdd1.weight
            temp_res=add(temp_bdd1,bdd2)
            the_successor.append([0,temp_res.weight,temp_res.node])
            the_successor+=[[succ[0],succ[1] * bdd1.weight,succ[2]] for succ in bdd1.node.successor[1:]]
        else:
            the_successor.append([0,bdd2.weight,bdd2.node])
            the_successor+=[[succ[0],succ[1]*bdd1.weight,succ[2]] for succ in bdd1.node.successor]
    else:
        the_key = k2
        if bdd2.node.successor[0][0]==0:
            temp_bdd2=BDD(bdd2.node.successor[0][2])
            temp_bdd2.weight=bdd2.node.successor[0][1]*bdd2.weight
            temp_res=add(bdd1,temp_bdd2)
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

def reduce_degree(bdd,min_degree):
    if len(min_degree)==0:
        return bdd.self_copy()
    
    max_index=max([global_index_order[k] for k in min_degree])
    
    k1=bdd.node.key
    
    if global_index_order[k1]>max_index:
        return bdd.self_copy()
    
    w=bdd.weight
    bdd.weight=1
    
    if find_computed_table(['r',bdd,min_degree]):
        res=find_computed_table(['r',bdd,min_degree])
        res.weight*=w
        bdd.weight=w
        return res
    
    the_successor=[]
    
    if k1 in min_degree:
        for succ in bdd.node.successor:
            temp_res=reduce_degree(BDD(succ[2]),min_degree)
            the_successor.append([succ[0]-min_degree[k1],succ[1]*temp_res.weight,temp_res.node])
        res=normalize(k1,the_successor)
    else:
        for succ in bdd.node.successor:
            temp_res=reduce_degree(BDD(succ[2]),min_degree)
            the_successor.append([succ[0],succ[1]*temp_res.weight,temp_res.node])
        res=normalize(k1,the_successor)
    insert_2_computed_table(['r',bdd,min_degree],res)
    res.weight*=w
    bdd.weight=w    
    return res


def normalize_2_fun(bdd1,bdd2):
    #bdd2/bdd1,the result is like [bdd1,1,bdd2/bdd1]
    if bdd1.weight==0:
        return [bdd2.self_copy(),bdd1,get_one_state()]
    
    if bdd2.weight==0:
        return [bdd1.self_copy(),get_one_state(),bdd2.self_copy()]    
    
    if bdd1.node.key==-1:
        temp=bdd2.self_copy()
        temp.weight/=bdd1.weight
        return [bdd1.self_copy(),get_one_state(),temp]
    
    if bdd1.node==bdd2.node:
        temp=get_one_state()
        temp.weight=bdd2.weight/bdd1.weight
        return [bdd1.self_copy(),get_one_state(),temp]
    
    # return [get_one_state(),bdd1,bdd2]

    w1=bdd1.weight
    w2=bdd2.weight
    bdd1.weight=1
    bdd2.weight=1
    
    if find_computed_table(['/',bdd1,bdd2]):
        [a,b,c]=find_computed_table(['/',bdd1,bdd2])
        a.weight*=w1
        c.weight*=w2/w1
        bdd1.weight=w1
        bdd2.weight=w2
        return [a,b,c]
    
    min_degree=dict()
    for key in bdd1.node.index_degree:
        if key in bdd2.node.index_degree:
            if min(bdd1.node.index_degree[key],bdd2.node.index_degree[key])>0:
                min_degree[key]=min(bdd1.node.index_degree[key],bdd2.node.index_degree[key])
                
    a=get_one_state()                
    if len(min_degree)>0:
#         print(min_degree)
        b=reduce_degree(bdd1,min_degree)
        c=reduce_degree(bdd2,min_degree)
        min_degree_order=[[global_index_order[k],k] for k in min_degree]
        while len(min_degree_order)>0:
            key=min_degree_order.pop()[1]
            a=normalize(key,[[min_degree[key],1,a.node]])        
    else:
        b=bdd1.self_copy()
        c=bdd2.self_copy()
    
    a.weight=b.weight
    b.weight=1
    c.weight/=a.weight
#     w=max(bdd1.weight,bdd2.weight)
#     if not get_int_key(bdd1.weight)==(0,0):
#         w=bdd1.weight
#     else:
#         w=bdd2.weight
#     a=BDD(Find_Or_Add_Unique_table(-1))
#     a.weight=w
#     b=bdd1.self_copy()
#     b.weight=bdd1.weight/w
#     c=bdd2.self_copy()
#     c.weight=bdd2.weight/w
    insert_2_computed_table(['/',bdd1,bdd2],[a,b,c])
    
    a.weight*=w1
    c.weight*=w2/w1
    bdd1.weight=w1
    bdd2.weight=w2
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

def conjugate(bdd):
    """"change the position of x and y, x<y in this bdd
    """
    v=bdd.node.key
    if v==-1:
        res=bdd.self_copy()
        res.weight=bdd.weight.conjugate()      
        return res
    the_successor=[]
    for k in bdd.node.successor:
        temp_res=conjugate(BDD(k[2]))
        the_successor.append([k[0],k[1].conjugate()*temp_res.weight,temp_res.node])
    res=normalize(v,the_successor)
    return res


def var_sort (var):
    global global_index_order, inverse_global_index_order
    # print('BDD 849', var, inverse_global_index_order, global_index_order)
    idx=[global_index_order[k] for k in var]
    idx.sort( )
    # print('BDD 853', idx)
    var_sort= [inverse_global_index_order[k] for k in idx]

    return var_sort

def cont(mode,bdd1,bdd2):
    #找出哪些要輸出(cont)/保留(out)
    var_out1=set(bdd1.key_2_index.values())
    var_out2=set(bdd2.key_2_index.values())
    print('BDD 877', 
          '\n   bdd:',bdd1, bdd2, 
          '\n   bdd node:', bdd1.node.key,bdd2.node.key, 
          '\n   key_2_index', bdd1.key_2_index, bdd2.key_2_index)
    var_out=list(var_out1.union(var_out2))
    var_out=var_sort(var_out) #Index 已含有比較大小的函數
    print('BDD 880',
          '\n   var_out1:', var_out1, 
          '\n   var_out2:', var_out2,
          '\n   var_out:', var_out)
    idx_2_key={-1:-1}
    key_2_idx={-1:-1}
    
    for new_key, idx in enumerate(var_out):
        if idx != -1:
            idx_2_key[idx]=new_key
            key_2_idx[new_key]=idx

    print('BDD 890',
          '\n   idx_2_key', idx_2_key,
          '\n   key_2_idx:',  key_2_idx)
    #找key的對應
    key_2_new_key=[[],[]]
    
    def set_new_key(key_2_index,key_2_new_key):
        vars = var_sort(key_2_index.values())
        for v in vars:
            assert v in idx_2_key, 'v=%s not found in idx_2_key'%v
            key_2_new_key.append(idx_2_key[v])
        print('BDD 896', key_2_index, key_2_new_key)
    
    set_new_key(bdd1.key_2_index, key_2_new_key[0])
    set_new_key(bdd2.key_2_index, key_2_new_key[1])
    
    if mode =='mul':
        bdd=mul2(bdd1,bdd2,key_2_new_key)
    if mode =='add':
        bdd=add2(bdd1,bdd2,key_2_new_key)
    # bdd.index_set=var_out
    bdd.key_2_index=key_2_idx

    return bdd


def mul2(bdd1,bdd2,key_2_new_key):
    """The contraction of two bdds, var_cont is in the form [[4,1],[3,2]]"""
    global global_index_order ,inverse_global_index_order

    k1 = bdd1.node.key
    k2 = bdd2.node.key
    w1 = bdd1.weight
    w2 = bdd2.weight
    print('BDD 928', k1, k2, key_2_new_key)
    
    weig = w1 * w2

    if k1==k2==-1:
        term = Find_Or_Add_Unique_table(-1)
        weig = 0 if get_int_key(weig) == (0,0) else weig
        return BDD(term, weight = weig)        

    if k1==-1:
        if w1==0:
            return BDD(bdd1.node,weight=0)

        key_2_index = {new_key:inverse_global_index_order[new_key] for new_key in key_2_new_key[1] if new_key != -1}
        key_2_index[-1] = -1

        return BDD(bdd2.node, weight = weig, key_2_index = key_2_index)
            
    if k2==-1:
        if w2==0:
            return BDD(bdd2.node ,weight=0)
 
        key_2_index = {new_key:inverse_global_index_order[new_key] for new_key in key_2_new_key[0] if new_key != -1}
        key_2_index[-1] = -1
        return BDD(bdd1.node, weight = weig, key_2_index = key_2_index)
    
    bdd1.weight=1
    bdd2.weight=1
    
    temp_key_2_new_key=[]
    temp_key_2_new_key.append(tuple([k for k in key_2_new_key[0][:(k1+1)]]))
    temp_key_2_new_key.append(tuple([k for k in key_2_new_key[1][:(k2+1)]]))

    temp_key_2_new_key=[]
    temp_key_2_new_key.append(tuple([k for k in key_2_new_key[0]]))
    temp_key_2_new_key.append(tuple([k for k in key_2_new_key[1]]))

    bdd=find_computed_table(['*',bdd1,bdd2,temp_key_2_new_key])
    # bdd=find_computed_table(['*',bdd1,bdd2,key_2_new_key])
    if bdd:
        bdd.weight=bdd.weight*weig
        bdd1.weight=w1
        bdd2.weight=w2        
        return bdd
    
    """BDD2.node 為什麼要這樣？"""
    bdd=BDD(Node(-1),weight=0.0)


    print('BDD 983',temp_key_2_new_key ,k1,k2)
    
  
    cont_order0 =  key_2_new_key[0][k1] 
    cont_order1 =  key_2_new_key[1][k2] 

    #A successor is in this form [degree, weight, node]
    if cont_order0 == cont_order1:
        for succ1 in bdd1.node.successor:
            for succ2 in bdd2.node.successor:
                temp_bdd1=BDD(succ1[2], weight=succ1[1], key_2_index= bdd1.get_key_2_index_subdict(succ1[2]))
                temp_bdd2=BDD(succ2[2], weight=succ2[1], key_2_index= bdd2.get_key_2_index_subdict(succ2[2]))
        
                # temp_res = mul2(temp_bdd1,temp_bdd2, key_2_new_key)
                temp_res = cont('mul',temp_bdd1 , temp_bdd2)

                print('BDD 999', key_2_new_key, k1)
                if not succ1[0]+succ2[0]==0: #Top node 同sin or cos
                    if succ1[0]+succ2[0]==2 and k1 % 2 == 1: #判斷是sin^2
                        temp_res1=normalize(k1 -1,[[2,-1,Find_Or_Add_Unique_table(-1)]])
                        # temp_res1=BDD(Node(k1 + 1), key_2_index = {k1:inverse_global_index_order[key_2_new_key[0][k1]], k1+1:inverse_global_index_order[key_2_new_key[0][k1+1]],-1:-1})
                        # temp_res1.node.successor = [[2,-1,Find_Or_Add_Unique_table(-1)]]
                        temp_res1.key_2_index = {k1:inverse_global_index_order[key_2_new_key[0][k1]], k1-1:inverse_global_index_order[key_2_new_key[0][k1-1]],-1:-1}
                        # temp_res1=mul2(temp_res1,temp_res,key_2_new_key)
                        # temp_res=add2(temp_res,temp_res1,key_2_new_key)
                        temp_res1 = cont('mul',temp_res1, temp_res)
                        temp_res = cont('add', temp_res, temp_res1)
                    else:
                        temp_res=normalize(k1 ,[[succ1[0]+succ2[0],temp_res.weight,temp_res.node]])
                        key_2_index = {k1:inverse_global_index_order[key_2_new_key[0][k1]], k1+1:inverse_global_index_order[key_2_new_key[0][k1+1]],-1:-1}
                        temp_res.key_2_index =  key_2_index
                # bdd=add2(bdd,temp_res,key_2_new_key)
                bdd = cont('add', bdd, temp_res)
    elif cont_order0 <= cont_order1:
        the_successor=[]
        for succ1 in bdd1.node.successor:
            temp_bdd1=BDD(succ1[2], weight=succ1[1], key_2_index= bdd1.get_key_2_index_subdict(succ1[2]))
            # temp_res=mul2(temp_bdd1, bdd2, key_2_new_key)
            temp_res = cont('mul',temp_bdd1, bdd2)
            the_successor.append([succ1[0],temp_res.weight,temp_res.node])
            bdd=normalize(key_2_new_key[0][k1],the_successor)
            key_2_index = {new_key:inverse_global_index_order[new_key] for new_key in key_2_new_key[1] if new_key != -1}
            key_2_index[-1] = -1
            bdd.key_2_index = key_2_index
        
    else:
        the_successor=[]
        for succ2 in bdd2.node.successor:
            temp_bdd2=BDD(succ2[2],weight=succ2[1], key_2_index= bdd1.get_key_2_index_subdict(succ2[2]))
            # temp_res=mul2(bdd1, temp_bdd2 ,key_2_new_key)
            temp_res = cont('mul',temp_bdd2, bdd1)
            the_successor.append([succ2[0],temp_res.weight,temp_res.node])
            bdd=normalize(key_2_new_key[1][k2],the_successor)     
            key_2_index = {new_key:inverse_global_index_order[new_key] for new_key in key_2_new_key[0] if new_key != -1}
            key_2_index[-1] = -1
            bdd.key_2_index = key_2_index   
            
    insert_2_computed_table(['*',bdd1,bdd2,temp_key_2_new_key],bdd)
    bdd.weight=bdd.weight*weig
    bdd1.weight=w1
    bdd2.weight=w2
    print('BDD 1031',bdd)
    return bdd      


def add2(bdd1,bdd2,key_2_new_key):
    """The apply function of two bdds. Mostly, it is used to do addition here."""
    global global_index_order    
    
    if bdd1.weight==0:
        return bdd2.self_copy()

    
    if bdd2.weight==0:
        return bdd1.self_copy()
    
    k1=bdd1.node.key
    k2=bdd2.node.key
    
    temp_key_2_new_key=[]
    temp_key_2_new_key.append(tuple([k for k in key_2_new_key[0]]))
    temp_key_2_new_key.append(tuple([k for k in key_2_new_key[1]]))

    if find_computed_table(['+',bdd1,bdd2,temp_key_2_new_key]):
        return find_computed_table(['+',bdd1,bdd2,temp_key_2_new_key])
    
    the_successor=[]
    cont_order0 =  key_2_new_key[0][k1] 
    cont_order1 =  key_2_new_key[1][k2] 
    cont_order0 = cont_order0 if cont_order0 != -1 else float('inf')
    cont_order1 = cont_order1 if cont_order1 != -1 else float('inf')
    print('BDD 1041', k1, k2, cont_order0, cont_order1, key_2_new_key)

    if cont_order0 == cont_order1:
        the_key = key_2_new_key[0][k1]
        degree1=[succ[0] for succ in bdd1.node.successor]
        degree2=[succ[0] for succ in bdd2.node.successor]
        the_successor+=[[succ[0],succ[1] * bdd1.weight,succ[2]] for succ in bdd1.node.successor if not succ[0] in degree2]
        the_successor+=[[succ[0],succ[1] * bdd2.weight,succ[2]] for succ in bdd2.node.successor if not succ[0] in degree1]
        com_degree=[d for d in degree1 if d in degree2]
        for d in com_degree:
            pos1=degree1.index(d)
            succ1 = bdd1.node.successor[pos1]
            temp_bdd1=BDD(succ1[2], weight= succ1[1] * bdd1.weight , key_2_index= bdd1.get_key_2_index_subdict(succ1[2]))
            
            pos2=degree2.index(d)
            succ2 = bdd2.node.successor[pos2]
            temp_bdd2=BDD(succ2[2], weight= succ2[1] * bdd2.weight , key_2_index= bdd2.get_key_2_index_subdict(succ2[2]))
          
            # temp_res=add2(temp_bdd1,temp_bdd2,key_2_new_key)
            temp_res = cont('add', temp_bdd1,temp_bdd2)
            the_successor.append([d,temp_res.weight,temp_res.node])
    elif cont_order0 <= cont_order1:
        the_key = key_2_new_key[0][k1]
        print('BDD 1061', bdd1.node.successor, the_key)
        # if len(bdd1.node.successor)!=0 and 
        if  bdd1.node.successor[0][0] == 0:
            succ1 = bdd1.node.successor[0]
            temp_bdd1=BDD(succ1[2], weight=succ1[1] * bdd1.weight , key_2_index= bdd1.get_key_2_index_subdict(succ1[2]))
            # temp_res=add2(temp_bdd1,bdd2,key_2_new_key)
            temp_res = cont('add', temp_bdd1,bdd2)
            the_successor.append([0,temp_res.weight,temp_res.node])
            the_successor+=[[succ[0],succ[1]*bdd1.weight,succ[2]] for succ in bdd1.node.successor[1:]]
        else:
            the_successor.append([0,bdd2.weight,bdd2.node])
            the_successor+=[[succ[0],succ[1]*bdd1.weight,succ[2]] for succ in bdd1.node.successor]
    else:
        the_key = key_2_new_key[1][k2]
        print('BDD 1114', bdd2.node.successor, the_key)
        # if len(bdd2.node.successor)!=0 and 
        if bdd2.node.successor[0][0]==0:
            succ2 = bdd2.node.successor[0]
            temp_bdd2=BDD(succ2[2], weight= succ2[1] *bdd2.weight , key_2_index= bdd2.get_key_2_index_subdict(succ2[2]))
 
            # temp_res=add2(bdd1,temp_bdd2,key_2_new_key)
            temp_res = cont('add', bdd1,temp_bdd2)
            the_successor.append([0,temp_res.weight,temp_res.node])
            the_successor+=[[succ[0],succ[1]*bdd2.weight,succ[2]] for succ in bdd2.node.successor[1:]]
        else:
            the_successor.append([0,bdd1.weight,bdd1.node])
            the_successor+=[[succ[0],succ[1]*bdd2.weight,succ[2]] for succ in bdd2.node.successor]          
            
    res = normalize(the_key,the_successor)
    j=0 if cont_order0 <= cont_order1 else 1
    key_2_index = {new_key:inverse_global_index_order[new_key] for new_key in key_2_new_key[j] if new_key != -1}
    key_2_index[-1] = -1
    res.key_2_index = key_2_index   
    insert_2_computed_table(['+',bdd1,bdd2,temp_key_2_new_key],res)
    print('BDD 1106',res)
#     print('-----------')
    return res


from functools import lru_cache, wraps, update_wrapper

def hashed_args(func):
    @wraps(func)
    def wrapper(*args):
        # print(*tuple(sorted(args)))
        return func(*tuple(sorted(args, key= hash)))
    
    update_wrapper(wrapper, func, ('cache_info', 'cache_clear','cache_parameters')) #繼承lru_cache的屬性
    return wrapper

@hashed_args #為了達成交換率，先排序輸入的兩個bdd
@lru_cache(maxsize=None)
def inner_mul(bdd1, bdd2):
    """The contraction of two bdds, var_cont is in the form [[4,1],[3,2]]"""

    '''original key'''
    k1=bdd1.node.key
    k2=bdd2.node.key
    w1=bdd1.weight
    w2=bdd2.weight
    print('BDD 1131 ',k1, k2 ,w1 ,w2)

    if k1==k2==-1:
        weig=bdd1.weight * bdd2.weight
        res=BDD(weig)        
        if get_int_key(weig)==(0,0):
            res.weight=0
            return res
        else:
            res.weight=weig
            return res        

    if k1==-1:
        if w1==0:
            bdd=BDD(bdd1.node)
            bdd.weight=0
            return bdd

        bdd=BDD(bdd2.node)
        bdd.weight = w1 * w2
        return bdd
            
    if k2==-1:
        if w2==0:
            bdd=BDD(bdd2.node)
            bdd.weight=0
            return bdd        

        bdd=BDD(bdd1.node)
        bdd.weight=w1 * w2
        return bdd
    
    bdd1.weight=1
    bdd2.weight=1
                
    bdd = BDD(bdd2.node)
    bdd.weight = 0
    # print ('BDD 1167', bdd.node.key)
    if k1 == k2:
        for succ1 in bdd1.node.successor:
            for succ2 in bdd2.node.successor:
                temp_bdd1=BDD(succ1[2])
                temp_bdd1.weight=succ1[1]
                temp_bdd2=BDD(succ2[2])
                temp_bdd2.weight=succ2[1]
                temp_res= mul3(temp_bdd1,temp_bdd2)
                if not succ1[0]+succ2[0]==0:
                    if succ1[0]+succ2[0]==2 and k1 % 2 == 0:
                        temp_res1 = normalize( k1 + 1,[[2,-1,Find_Or_Add_Unique_table(-1)]])
                        temp_res1 = mul3(temp_res1,temp_res)
                        temp_res = add3(temp_res,temp_res1)
                    else:
                        temp_res = normalize( k1, [[succ1[0]+succ2[0],temp_res.weight,temp_res.node]])
                bdd=add3(bdd,temp_res)
                
    elif k1 <= k2:
        the_successor=[]
        for succ1 in bdd1.node.successor:
            temp_bdd1=BDD(succ1[2])
            temp_bdd1.weight=succ1[1]
            temp_res=mul3(temp_bdd1,bdd2)
            the_successor.append([succ1[0],temp_res.weight,temp_res.node])
            bdd=normalize( k1, the_successor)
    else:
        the_successor=[]
        for succ2 in bdd2.node.successor:
            temp_bdd2=BDD(succ2[2])
            temp_bdd2.weight=succ2[1]
            temp_res=mul3(bdd1,temp_bdd2)
            the_successor.append([succ2[0],temp_res.weight,temp_res.node])
            bdd=normalize( k2,the_successor)        
            
    bdd.weight=bdd.weight * w1 * w2
    bdd1.weight=w1
    bdd2.weight=w2
    return bdd      


@hashed_args #為了達成交換率，先排序輸入的兩個bdd
@lru_cache(maxsize=None)
def mul3(bdd1, bdd2):
    print('BDD 1202',bdd1.key_2_index)
    print('BDD 1203',bdd2.key_2_index)

    #get union key_2_index and index_set
    index_set=list(set(bdd1.key_2_index.values()).union(set(bdd2.key_2_index.values())))
    index_set.sort
    key_2_index=dict()
    for k in index_set:
        key_2_index[index_set.index(k)] = inverse_global_index_order[k]
    key_2_index[-1] = -1

    #How to change bdd and input to inner mul? 
    
    print('BDD 1217',bdd1.node.key, bdd2.node.key)
    print('BDD 1218',bdd1.key_2_index[bdd1.node.key], bdd2.key_2_index[bdd2.node.key])
    
    bdd=inner_mul(bdd1, bdd2)

    #change back the idex_set
    # bdd.node.key = 
    bdd.index_set = index_set
    bdd.key_2_index = key_2_index

    print('BDD 1208', index_set, key_2_index)
    return bdd

@hashed_args #為了達成交換率，先排序輸入的兩個bdd
@lru_cache(maxsize=None)
def add3(bdd1, bdd2):
    """The apply function of two bdds. Mostly, it is used to do addition here."""
    global global_index_order    

    if bdd1.weight == 0:
        return bdd2.self_copy()
    
    if bdd2.weight == 0:
        return bdd1.self_copy()
    
    if bdd1.node == bdd2.node:
        weig = bdd1.weight + bdd2.weight
        if get_int_key(weig) == (0,0):
            term = Find_Or_Add_Unique_table(-1)
            bdd = BDD(term)
            bdd.weight=0
            return bdd
        else:
            bdd = BDD(bdd1.node)
            bdd.weight = weig
            return bdd
    
    k1 = bdd1.node.key
    k2 = bdd2.node.key
    print('BDD 1218', k1, k2)
    the_successor = []
    if k1 == k2:
        the_key = k1
        degree1 = [succ[0] for succ in bdd1.node.successor]
        degree2 = [succ[0] for succ in bdd2.node.successor]
        the_successor += [[succ[0],succ[1] * bdd1.weight,succ[2]] for succ in bdd1.node.successor if not succ[0] in degree2]
        the_successor += [[succ[0],succ[1] * bdd2.weight,succ[2]] for succ in bdd2.node.successor if not succ[0] in degree1]
        com_degree=[d for d in degree1 if d in degree2]
        for d in com_degree:
            pos1=degree1.index(d)
            temp_bdd1=BDD(bdd1.node.successor[pos1][2])
            temp_bdd1.weight=bdd1.node.successor[pos1][1]*bdd1.weight
            pos2=degree2.index(d)
            temp_bdd2=BDD(bdd2.node.successor[pos2][2])
            temp_bdd2.weight=bdd2.node.successor[pos2][1]*bdd2.weight          
            temp_res=add3(temp_bdd1,temp_bdd2)
#             print(temp_bdd1,temp_bdd2,temp_res)
            the_successor.append([d,temp_res.weight,temp_res.node])
    elif k1 <= k2:
    # elif global_index_order[k1]<=global_index_order[k2]:
        the_key = k1
        print('BDD 649',bdd1.node)
        if bdd1.node.successor[0][0] == 0:
            temp_bdd1=BDD(bdd1.node.successor[0][2])
            temp_bdd1.weight=bdd1.node.successor[0][1] * bdd1.weight
            temp_res=add3(temp_bdd1,bdd2)
            the_successor.append([0,temp_res.weight,temp_res.node])
            the_successor+=[[succ[0],succ[1] * bdd1.weight,succ[2]] for succ in bdd1.node.successor[1:]]
        else:
            the_successor.append([0,bdd2.weight,bdd2.node])
            the_successor+=[[succ[0],succ[1] * bdd1.weight,succ[2]] for succ in bdd1.node.successor]
    else:
        the_key = k2
        if bdd2.node.successor[0][0]==0:
            temp_bdd2=BDD(bdd2.node.successor[0][2])
            temp_bdd2.weight=bdd2.node.successor[0][1] * bdd2.weight
            temp_res=add3(bdd1,temp_bdd2)
            the_successor.append([0,temp_res.weight,temp_res.node])
            the_successor+=[[succ[0],succ[1] * bdd2.weight,succ[2]] for succ in bdd2.node.successor[1:]]
        else:
            the_successor.append([0,bdd1.weight,bdd1.node])
            the_successor+=[[succ[0],succ[1] * bdd2.weight,succ[2]] for succ in bdd2.node.successor]          
            
    bdd = normalize(the_key,the_successor)
    insert_2_computed_table(['+',bdd1,bdd2],bdd)
#     print(bdd)
#     print('-----------')
    return bdd