import numpy as np
import copy
import time
import random
from graphviz import Digraph
from IPython.display import Image
from sympy import *

"""Define global variables"""
computed_table = dict()
unique_table = dict()
global_index_order = dict()
global_node_idx=0
add_find_time=0
add_hit_time=0
cont_find_time=0
cont_hit_time=0
epi=0.000001

    
class Node:
    """To define the node of TDD"""
    def __init__(self,key,num=2):
        self.idx = 0
        self.key = key
        self.succ_num=2
        self.out_weight=[1]*num
        self.successor=[None]*num
        self.expr = None

class BDD:
    def __init__(self,node):
        """BDD"""
        self.weight=1

        if isinstance(node,Node):
            self.node=node
        else:
            self.node=Node(node)
            
            
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
        dot.node('-0','',shape='none')
        dot.edge('-0',str(self.node.idx),color="blue",label=str(complex(round(self.weight.real,2),round(self.weight.imag,2))))
        dot.format = 'png'
        return Image(dot.render('output'))
        
    def __eq__(self,other):
        if self.node==other.node and get_int_key(self.weight)==get_int_key(other.weight):
            return True
        else:
            return False
        
    def __add__(self, g):        
        return add(self,g)
    
    def __mul__(self, g):
        return mul(self,g)
    
    def self_normalize(self, g):
        return normalize_2_fun(g,self)
    
    def expr(self):
        if self.node.key==-1:
            return self.weight
        res=get_expr(self.node)
        return self.weight*res
    def __repr__(self):
        return str(self.expr())
    
def layout(node,dot=Digraph(),succ=[]):
    col=['red','blue','black','green']
    dot.node(str(node.idx), str(node.key), fontname="helvetica",shape="circle",color="red")
    for k in range(node.succ_num):
        if node.successor[k]:
            label1=str(complex(round(node.out_weight[k].real,2),round(node.out_weight[k].imag,2)))
            if not node.successor[k] in succ:
                dot=layout(node.successor[k],dot,succ)
                dot.edge(str(node.idx),str(node.successor[k].idx),color=col[k%4],label=label1)
                succ.append(node.successor[k])
            else:
                dot.edge(str(node.idx),str(node.successor[k].idx),color=col[k%4],label=label1)
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
    if not node in node_set:
        node_set.add(node)
        for k in range(node.succ_num):
            if node.successor[k]:
                node_set = get_node_set(node.successor[k],node_set)
    return node_set

def Find_Or_Add_Unique_table(x,weigs=[],succ_nodes=[]):
    """To return a node if it already exist, creates a new node otherwise"""
    global global_node_idx,unique_table
    
    if x==-1:
        if unique_table.__contains__(x):
            return unique_table[x]
        else:
            res=Node(x)
            res.idx=0
            unique_table[x]=res
#             print('bbbbb:',res)
        return res
    temp_key=[x]
    for k in range(len(weigs)):
        temp_key.append(get_int_key(weigs[k]))
        temp_key.append(succ_nodes[k])
    temp_key=tuple(temp_key)
    if temp_key in unique_table:
        return unique_table[temp_key]
    else:
        res=Node(x,len(succ_nodes))
        global_node_idx+=1
        res.idx=global_node_idx
        res.out_weight=weigs
        res.successor=succ_nodes
        unique_table[temp_key]=res
    return res


def normalize(x,the_successors):
    """The normalize and reduce procedure"""
    global epi
    all_equal=True
    for k in range(1,len(the_successors)):
        if the_successors[k]!=the_successors[0]:
            all_equal=False
            break
    if all_equal:
        return the_successors[0]
    
    weigs=[succ.weight for succ in the_successors]
    
    weigs_abs=[np.around(abs(weig)/epi) for weig in weigs]
    weig_max=weigs[weigs_abs.index(max(weigs_abs))]
    weigs=[weig/weig_max for weig in weigs]
    for k in range(len(the_successors)):
        if get_int_key(weigs[k])==(0,0):
            node=Find_Or_Add_Unique_table(-1)
            the_successors[k]=BDD(node)
            the_successors[k].weight=weigs[k]=0    
    succ_nodes=[succ.node for succ in the_successors]
    node=Find_Or_Add_Unique_table(x,weigs,succ_nodes)
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
    else:
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
    else:
        the_key = ('/',get_int_key(item[1].weight),item[1].node,get_int_key(item[2].weight),item[2].node)
        computed_table[the_key] = (res[0].weight,res[0].node,res[1].weight,res[1].node,res[2].weight,res[2].node)
    

def get_bdd(f):
    global global_index_order 
    if isinstance(f,int) or isinstance(f,float) or isinstance(f,complex):
        tdd=get_one_state()
        tdd.weight=f
        return tdd
    
    try:
        sm=[str(item) for item in f.free_symbols]
        res=[]
        for k in sm:
            if k[1]=='n':
                t=k[0]+k[2:]
            else:
                t=k
            tt=(global_index_order[t],t)
            if not tt in res:
                res.append(tt)
        res.sort()
    except:
        tdd=get_one_state()
        tdd.weight=complex(f)
        return tdd
    
    key=res[0][1]
    key2=res[0][1][0]+'n'+res[0][1][1:]
    f1=f.xreplace({symbols(key):0,symbols(key2):1})
    f2=f.xreplace({symbols(key):1,symbols(key2):0})
#     print(sm,f1,f2)
    bdd1=get_bdd(f1)
    bdd2=get_bdd(f2)
    bdd=normalize(key,[bdd1,bdd2])
    return bdd


def mul(tdd1,tdd2):
    """The contraction of two TDDs, var_cont is in the form [[4,1],[3,2]]"""
    global global_index_order 
    
    k1=tdd1.node.key
    k2=tdd2.node.key
    w1=tdd1.weight
    w2=tdd2.weight
    
    if tdd1.node==tdd2.node:
        weig=tdd1.weight*tdd2.weight
        if get_int_key(weig)==(0,0):
            term=Find_Or_Add_Unique_table(-1)
            res=BDD(term)
            res.weight=0
            return res
        else:
            res=BDD(tdd1.node)
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
                
    if global_index_order[k1]<=global_index_order[k2]:
        the_key=k1
    else:
        the_key=k2
        
        
    the_successors=[]
    for k in range(tdd1.node.succ_num):
        res=mul(Slicing(tdd1,the_key,k),Slicing(tdd2,the_key,k))
        the_successors.append(res)
            
    tdd=normalize(the_key,the_successors)
    insert_2_computed_table(['*',tdd1,tdd2],tdd)
    tdd.weight=tdd.weight*w1*w2

    tdd1.weight=w1
    tdd2.weight=w2
    
    return tdd
    
def Slicing(tdd,x,c):
    """Slice a TDD with respect to x=c"""
    global global_index_order 
    k=tdd.node.key
    
    if k==-1:
        return tdd.self_copy()
    
    if global_index_order[k]>global_index_order[x]:
        return tdd.self_copy()
    
    if k==x:
        res=BDD(tdd.node.successor[c])
        res.weight=tdd.node.out_weight[c]
        return res
    else:
        print("Not supported yet!!!")
        
        
def Slicing2(tdd,x,c):
    """Slice a TDD with respect to x=c"""
    global global_index_order 
    k=tdd.node.key
    
    if k==-1:
        return tdd.self_copy()
    
    if global_index_order[k]>global_index_order[x]:
        return tdd.self_copy()
    
    if k==x:
        res=BDD(tdd.node.successor[c])
        res.weight=tdd.node.out_weight[c]*tdd.weight
        return res
    else:
        print("Not supported yet!!!")        
        

def add(tdd1,tdd2):
    """The apply function of two TDDs. Mostly, it is used to do addition here."""
    global global_index_order    

    k1=tdd1.node.key
    k2=tdd2.node.key
#     print(tdd1,k1)
#     print(tdd2,k2)
#     print(tdd1.node.key,tdd2.node.key,tdd1.node,tdd2.node,tdd1.node.out_weight,tdd2.node.out_weight,tdd1.node.successor,tdd2.node.successor)
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
    
    if global_index_order[k1]<=global_index_order[k2]:
        the_key=k1
    else:
        the_key=k2
        
    the_successors=[]
    

    for k in range(2):
        res=add(Slicing2(tdd1,the_key,k),Slicing2(tdd2,the_key,k))
        the_successors.append(res)
            
    res = normalize(the_key,the_successors)
    insert_2_computed_table(['+',tdd1,tdd2],res)
    return res


def normalize_2_fun(tdd1,tdd2):
    #tdd2/tdd1
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
    if 0:
#     if tdd1.node.out_weight[0]!=0 and tdd1.node.out_weight[1]!=0:
        print(tdd1.node.out_weight)
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
#             print('11111')
            return [a,b,c]

        k1=tdd1.node.key
        k2=tdd2.node.key
    
        if global_index_order[k1]<=global_index_order[k2]:
            the_key=k1
        else:
            the_key=k2
        
        [a1,b1,c1]=normalize_2_fun(Slicing(tdd1,the_key,0),Slicing(tdd2,the_key,0))
        [a2,b2,c2]=normalize_2_fun(Slicing(tdd1,the_key,1),Slicing(tdd2,the_key,1))
    
        a= normalize(the_key,[a1,a2])
        b= normalize(the_key,[b1,b2])
        c= normalize(the_key,[c1,c2])
        insert_2_computed_table(['/',tdd1,tdd2],[a,b,c])

        a.weight*=w1
        c.weight*=w2/w1    
        tdd1.weight=w1
        tdd2.weight=w2
        return [a,b,c]
    else:
        if find_computed_table(['/',tdd1,tdd2]):
            [a,b,c]=find_computed_table(['/',tdd1,tdd2])
#             print('11111')
            return [a,b,c]

        k1=tdd1.node.key
        k2=tdd2.node.key
    
        if global_index_order[k1]<=global_index_order[k2]:
            the_key=k1
        else:
            the_key=k2
        
        [a1,b1,c1]=normalize_2_fun(Slicing2(tdd1,the_key,0),Slicing2(tdd2,the_key,0))
        [a2,b2,c2]=normalize_2_fun(Slicing2(tdd1,the_key,1),Slicing2(tdd2,the_key,1))
    
        a= normalize(the_key,[a1,a2])
        b= normalize(the_key,[b1,b2])
        c= normalize(the_key,[c1,c2])
        insert_2_computed_table(['/',tdd1,tdd2],[a,b,c])
        return [a,b,c]        
    
    
def get_expr(node):
    if node.key==-1:
        return 1
    if node.expr:
        return node.expr
    x=node.key
    xn=x[0]+'n'+x[1:]
    l=get_expr(node.successor[0])
    r=get_expr(node.successor[1])
    if l==r:
        res=(node.out_weight[1]*symbols(x)+node.out_weight[0]*symbols(xn))*l
    else:
        res=node.out_weight[1]*symbols(x)*l+node.out_weight[0]*symbols(xn)*r
    node.expr=res
    return res