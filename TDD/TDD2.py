import numpy as np
import copy
import time
import random
from sympy import *
from sympy.parsing.sympy_parser import parse_expr
from graphviz import Digraph
from IPython.display import Image

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

class Index:
    """The index, here idx is used when there is a hyperedge"""
    def __init__(self,key,idx=0):
        self.key = key
        self.idx = idx
        
    def __eq__(self,other):
        if self.key == other.key and self.idx == other.idx:
            return True
        else:
            return False
        
    def __lt__(self,other):
        if global_index_order[self.key] < global_index_order[other.key]:
            return True
        elif self.key == other.key and self.idx<other.idx:
            return True
        
        return False
    
    def __str__(self):
        return str((self.key,self.idx))

    
class Node:
    """To define the node of TDD"""
    def __init__(self,key,num=2):
        self.idx = 0
        self.key = key
        self.succ_num=num
        self.out_weight=[1]*num
        self.successor=[None]*num
        self.meas_prob=[]


class TDD:
    def __init__(self,node):
        """TDD"""
        self.weight=[1]
        
        self.index_set=[]
        
        self.key_2_index=dict()
        self.index_2_key=dict()
        self.key_width=dict() #Only used when change TDD to numpy
        
        if isinstance(node,Node):
            self.node=node
        else:
            self.node=Node(node)
            
            
    def node_number(self):
        node_set=set()
        node_set=get_node_set(self.node,node_set)
        return len(node_set)
    
    def self_copy(self):
        temp = TDD(self.node)
        temp.weight = self.weight
        temp.index_set = copy.copy(self.index_set)
        temp.key_2_index=copy.copy(self.key_2_index)
        temp.index_2_key=copy.copy(self.index_2_key)
        return temp
    
    def show(self,real_label=True):
        edge=[]              
        dot=Digraph(name='reduced_tree')
        dot=layout(self.node,self.key_2_index,dot,edge,real_label)
        dot.node('-0','',shape='none')
        dot.edge('-0',str(self.node.idx),color="blue",label=str(self.weight))
        dot.format = 'png'
        return Image(dot.render('output'))
        
    def __eq__(self,other):
        if self.node==other.node and self.weight==other.weight:
            return True
        else:
            return False
        
def layout(node,key_2_idx,dot=Digraph(),succ=[],real_label=True):
    col=['red','blue','black','green']
    if real_label and node.key in key_2_idx:
        if node.key==-1:
            dot.node(str(node.idx), str(1), fontname="helvetica",shape="circle",color="red")
        else:
            dot.node(str(node.idx), key_2_idx[node.key], fontname="helvetica",shape="circle",color="red")
    else:
        dot.node(str(node.idx), str(node.key), fontname="helvetica",shape="circle",color="red")
    for k in range(node.succ_num):
        if node.successor[k]:
            label1=str(node.out_weight[k])
            if not node.successor[k] in succ:
                dot=layout(node.successor[k],key_2_idx,dot,succ,real_label)
                dot.edge(str(node.idx),str(node.successor[k].idx),color=col[k%4],label=label1)
                succ.append(node.successor[k])
            else:
                dot.edge(str(node.idx),str(node.successor[k].idx),color=col[k%4],label=label1)
    return dot        


def to_cnf_one_term(expr,s):
    sm=[str(item) for item in expr.free_symbols]
#     print(sm)
    sl=str(expr).split('+')
#     print(sl)
    l=[]
    r=[]
    if s in sm:
        for item in sl:
            if s in item:
                l.append(parse_expr(item)/parse_expr(s))
            else:
                r.append(parse_expr(item))
    sl0=str(l[0]).split('*')
    c=1
    for item in sl0:
        if not item in sm:
            try:
                c*=parse_expr(item)
            except:
                pass
#     print(c)
    try:
#         if get_int_key(float(c))==(0,0):
#             c=0
#             l=[item*0 for item in l]
#         else:
            l=[item/c for item in l]
    except:
        l=[item/c for item in l]
#     print(type(c))
#     print(c,l)
    return c,l
            


def to_cnf2(expr,n=5):
    expr=nsimplify(expr,tolerance=1e-6,rational=False)
    res=[]
    temp=factor_list(expr)
#     print('acc:',temp)
    for item in temp[1]:
        if item:
            res.append(item[0])
    if res:
        res[0]*=temp[0]
    else:
        res=[simplify(0*symbols('x'))]
#     print(res)
    
    return res

# def to_cnf2(expr,n=5):
#     sm=[str(item) for item in expr.free_symbols]
# #     print(sm)
#     sl=str(expr).split('+')
# #     print(sl)
#     if len(sl)<=1:
#         return [expr]
    
#     for k in range(n):
#         x='x'+str(k+1)
#         xn='xn'+str(k+1)
#         if not x in sm and not xn in sm:
#             continue
#         if x in sm and not xn in sm:
#             all_have_x = True
#             for item in sl:
#                 if not x in item:
#                     all_have_x=False
#                     break
#             if all_have_x:
#                 c1,l1 = to_cnf_one_term(expr,x)
#                 new_expr=0
#                 for item in l1:
#                     new_expr+=item
#                 r=to_cnf2(new_expr,n)
#                 r.insert(0,c1*parse_expr(x))
#                 return r
#             else:
#                 continue
#         if xn in sm and not x in sm:
#             all_have_xn = True
#             for item in sl:
#                 if not xn in item:
#                     all_have_xn=False
#                     break
#             if all_have_xn:
#                 c1,l1 = to_cnf_one_term(expr,xn)
#                 new_expr=0
#                 for item in l1:
#                     new_expr+=item
#                 r=to_cnf2(new_expr,n)
#                 r.insert(0,c1*parse_expr(xn))
#                 return r
#             else:
#                 continue
#         all_have_x_or_xn = True
#         for item in sl:
#             if not xn in item and not x in item:
#                 all_have_x_or_xn=False
#                 break
#         if all_have_x_or_xn:
#             c1,l1 = to_cnf_one_term(expr,x)
#             c2,l2 = to_cnf_one_term(expr,xn)
#             if not l1==l2:
#                 continue
#             else:
#                 new_expr=0
#                 for item in l1:
#                     new_expr+=item
#                 r=to_cnf2(new_expr,n)
#                 r.insert(0,c1*parse_expr(x)+c2*parse_expr(xn))
#                 return r
#         else:
#             continue
#     return [expr]



        
def Ini_TDD(index_order=[]):
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
    return get_identity_tdd()

def Clear_TDD():
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
#     return 1

def get_identity_tdd():
    node = Find_Or_Add_Unique_table(-1)
    tdd = TDD(node)
    tdd.index_2_key={-1:-1}
    tdd.key_2_index={-1:-1}
    return tdd

def get_unique_table():
    return global_unique_table

def get_unique_table_num():
    return len(global_unique_table)

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

def get_weight(sl):
    if not isinstance(sl,list):
        return sl
    r=1
    for item in sl:
        r*=item
    return simplify(r)

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
        return res
    temp_key=[x]
    for k in range(len(weigs)):
        temp_key.append(get_weight(weigs[k]))
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


def merge_I(expr):
    sm=[str(item) for item in expr.free_symbols]
    if len(sm)==0:
        return expr
    
    sl=str(expr).split('+')
    if len(sl)<=1:
        return expr
    res=0
    
    for s in sm:
        w=0
        for item in sl:
            if s in item:
                w+=parse_expr(item)/parse_expr(s)
        res+=w*parse_expr(s)
    return res

def get_a_weight(expr,s):
    sm=[str(item) for item in expr.free_symbols]
    if not s in sm:
        return 0
    sl=str(expr).split('+')
    w=0
    for item in sl:
        if s in item:
            w+=simplify(parse_expr(item)/parse_expr(s))
    print('wwwwwwwww\n',simplify(w))
    return simplify(w)

def div_expr(expr_t1,expr_t2):
    sm=[str(item) for item in expr_t1.free_symbols]
    sm.sort()    
    w0=get_a_weight(expr_t1,sm[0])
    w1=get_a_weight(expr_t1,sm[1])
    w20=get_a_weight(expr_t2,sm[0])
    w21=get_a_weight(expr_t2,sm[1])
    return simplify(w0/w20)*parse_expr(sm[0])+simplify(w1/w21)*parse_expr(sm[1])


def norm2(expr1,expr2):
#     print('ccc:',expr1,expr2)
    res=[]
    for expr_t1 in expr2:
        sm=[str(item) for item in expr_t1.free_symbols]
        sm.sort()
        if sm:
            for expr_t2 in expr1:
                sm2=[str(item) for item in expr_t2.free_symbols]
                sm2.sort()
                if sm==sm2:
                    expr=div_expr(expr_t1,expr_t2)
                    res.append(expr)
        else:
            res.append(expr_t1)
    if not res:
        return expr1
    return res


def renorm(tdd):
    if tdd.node.key==-1:
        return tdd
    
    l=renorm(TDD(tdd.node.successor[0]))
    r=renorm(TDD(tdd.node.successor[1]))
    if l.weight==[1]:
        l.weight=tdd.node.out_weight[0]
    else:
        l.weight+=tdd.node.out_weight[0]
    if r.weight==[1]:        
        r.weight=tdd.node.out_weight[1]
    else:
        r.weight+=tdd.node.out_weight[1]        
        
    the_successors=[l,r]
    
    return normalize2(tdd.node.key,the_successors)

def normalize2(x,the_successors):
    """The normalize and reduce procedure"""
    
    weigs=[succ.weight for succ in the_successors]
    print('w:',weigs)
    weig=[]
    for item in weigs[0]:
        if item in weigs[1]:
            weig.append(item)
            weigs[0].remove(item)
            weigs[1].remove(item)
            
    if weigs[0] and weigs[1]:
        for k in weigs[0]:
            weig.append(k)
        
        weigs2=norm2(weigs[0],weigs[1])
        weigs=[[simplify(symbols('x')/symbols('x'))],weigs2]
    print('weigs2:',weigs)
    succ_nodes=[succ.node for succ in the_successors]
    node=Find_Or_Add_Unique_table(x,weigs,succ_nodes)
    res=TDD(node)
    if len(weig)==0:
        res.weight=[1]
    else:
        res.weight=weig
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
#     print('bbb:',weigs)
    weig=[]
    for item in weigs[0]:
        if item in weigs[1]:
            weig.append(item)
            weigs[0].remove(item)
            weigs[1].remove(item)
            
#     if weigs[0] and weigs[1]:
#         for k in weigs[0]:
#             weig.append(k)
        
#         weigs2=norm2(weigs[0],weigs[1])
#         weigs=[[simplify(symbols('x')/symbols('x'))],weigs2]
        
    succ_nodes=[succ.node for succ in the_successors]
    node=Find_Or_Add_Unique_table(x,weigs,succ_nodes)
    res=TDD(node)
    if len(weig)==0:
        res.weight=[1]
    else:
        res.weight=weig
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
        the_key=('s',item[1].weight,item[1].node,temp_key,item[3])
        if computed_table.__contains__(the_key):
            res = computed_table[the_key]
            tdd = TDD(res[1])
            tdd.weight = res[0]
            return tdd
    elif item[0] == '+':
        the_key=('+',get_weight(item[1].weight),item[1].node,get_weight(item[2].weight),item[2].node)
        add_find_time+=1
        if computed_table.__contains__(the_key):
            res = computed_table[the_key]
            tdd = TDD(res[1])
            tdd.weight = res[0]
            add_hit_time+=1
            return tdd
        the_key=('+',get_weight(item[2].weight),item[2].node,get_weight(item[1].weight),item[1].node)
        if computed_table.__contains__(the_key):
            res = computed_table[the_key]
            tdd = TDD(res[1])
            tdd.weight = res[0]
            add_hit_time+=1
            return tdd
    else:
        the_key=('*',get_weight(item[1].weight),item[1].node,get_weight(item[2].weight),item[2].node,item[3][0],item[3][1],item[4])
        cont_find_time+=1
        if computed_table.__contains__(the_key):
            res = computed_table[the_key]
            tdd = TDD(res[1])
            tdd.weight = res[0]
            cont_hit_time+=1            
            return tdd
        the_key=('*',get_weight(item[2].weight),item[2].node,get_weight(item[1].weight),item[1].node,item[3][1],item[3][0],item[4])
        if computed_table.__contains__(the_key):
            res = computed_table[the_key]
            tdd = TDD(res[1])
            tdd.weight = res[0]
            cont_hit_time+=1            
            return tdd
    return None

def insert_2_computed_table(item,res):
    """To insert an item to the computed table"""
    global computed_table,cont_time,find_time,hit_time
    if item[0]=='s':
        temp_key=item[1].index_2_key[item[2]]
        the_key = ('s',item[1].weight,item[1].node,temp_key,item[3])
    elif item[0] == '+':
        the_key = ('+',get_weight(item[1].weight),item[1].node,get_weight(item[2].weight),item[2].node)
    else:
        the_key = ('*',get_weight(item[1].weight),item[1].node,get_weight(item[2].weight),item[2].node,item[3][0],item[3][1],item[4])
    computed_table[the_key] = (res.weight,res.node)
    
def get_index_2_key(var):
    var_sort=copy.copy(var)
    var_sort.sort()
    var_sort.reverse()
    idx_2_key={-1:-1}
    key_2_idx={-1:-1}
    n=0
    for idx in var_sort:
        if not idx.key in idx_2_key:
            idx_2_key[idx.key]=n
            key_2_idx[n]=idx.key
            n+=1
    return idx_2_key,key_2_idx
    
def get_tdd(U,var=[]):
    
#     if len(var)==0 or not isinstance(var[0],Index):
#         return np_2_tdd(U,var)
    
    idx_2_key, key_2_idx = get_index_2_key(var)
    order=[]
    for idx in var:
        order.append(idx_2_key[idx.key])
    tdd = np_2_tdd(U,order)
    tdd.index_2_key=idx_2_key
    tdd.key_2_index=key_2_idx
    tdd.index_set=var
    
#     if not order:
#         order=list(range(U_dim))
        
#     for k in range(max(order)+1):
#         split_pos=order.index(k)
#         tdd.key_width[k]=U.shape[split_pos]
        
    return tdd    
    

def np_2_tdd(U,order=[],key_width=True):
    #index is the index_set as the axis order of the matrix
    U_dim=U.ndim
    U_shape=U.shape
    if sum(U_shape)==U_dim:
        node=Find_Or_Add_Unique_table(-1)
        res=TDD(node)
        for k in range(U_dim):
            U=U[0]
        res.weight=[simplify(U)]
        return res
    
    if not order:
        order=list(range(U_dim))
    
    if key_width:
        the_width=dict()
        for k in range(max(order)+1):
            split_pos=order.index(k)
            the_width[k]=U.shape[split_pos]
            
    x=max(order)
    split_pos=order.index(x)
    order[split_pos]=-1
    split_U=np.split(U,U_shape[split_pos],split_pos)

    while x in order:
        split_pos=order.index(x)
        for k in range(len(split_U)):
            split_U[k]=np.split(split_U[k],U_shape[split_pos],split_pos)[k]
        order[split_pos]=-1
    
    the_successors=[]
    for k in range(U_shape[split_pos]):
        res=np_2_tdd(split_U[k],copy.copy(order),False)
        the_successors.append(res)
    tdd = normalize(x,the_successors)
    
    if key_width:
        tdd.key_width=the_width
    
    return tdd


    
def cont(tdd1,tdd2):

    var_cont=[var for var in tdd1.index_set if var in tdd2.index_set]
    var_out1=[var for var in tdd1.index_set if not var in var_cont]
    var_out2=[var for var in tdd2.index_set if not var in var_cont]

    var_out=var_out1+var_out2
    var_out.sort()
    var_out_idx=[var.key for var in var_out]
    var_cont_idx=[var.key for var in var_cont]
    var_cont_idx=[var for var in var_cont_idx if not var in var_out_idx]
    
    idx_2_key={-1:-1}
    key_2_idx={-1:-1}
    
    n=0
    for k in range(len(var_out_idx)-1,-1,-1):
        if not var_out_idx[k] in idx_2_key:
            idx_2_key[var_out_idx[k]]=n
            key_2_idx[n]=var_out_idx[k]
            n+=1
        
    key_2_new_key=[[],[]]
    cont_order=[[],[]]
    for k in range(len(tdd1.key_2_index)-1):
        v=tdd1.key_2_index[k]
        if v in idx_2_key:
            key_2_new_key[0].append(idx_2_key[v])
        else:
            key_2_new_key[0].append('c')
        cont_order[0].append(global_index_order[v])
        
    cont_order[0].append(float('inf'))
    
    for k in range(len(tdd2.key_2_index)-1):     
        v=tdd2.key_2_index[k]
        if v in idx_2_key:
            key_2_new_key[1].append(idx_2_key[v])
        else:
            key_2_new_key[1].append('c')
        cont_order[1].append(global_index_order[v])
    cont_order[1].append(float('inf'))

    tdd=contract(tdd1,tdd2,key_2_new_key,cont_order,len(set(var_cont_idx)))
    tdd.index_set=var_out
    tdd.index_2_key=idx_2_key
    tdd.key_2_index=key_2_idx
    key_width=dict()
    for k1 in range(len(key_2_new_key[0])):
        if not key_2_new_key[0][k1]=='c' and not key_2_new_key[0][k1] ==-1:
            key_width[key_2_new_key[0][k1]]=tdd1.key_width[k1]
    for k2 in range(len(key_2_new_key[1])):
        if not key_2_new_key[1][k2]=='c' and not key_2_new_key[1][k2] ==-1:
            key_width[key_2_new_key[1][k2]]=tdd2.key_width[k2]             
   
    tdd.key_width=key_width
#     print(tdd1.key_width,tdd2.key_width,tdd.key_width)
    return tdd


def mul_weight(w1,w2):
    if w1==[0] or w2==[0]:
        return [0]
    if w1==[1]:
        return w2
    if w2==[1]:
        return w1
    return to_cnf2(get_sum_form(w1+w2))

# def mul_weight(w1,w2):
#     if w1==[0] or w2==[0]:
#         return [0]
#     if w1==[1]:
#         return w2
#     if w2==[1]:
#         return w1
#     return w1+w2

def get_sum_form(sl):
#     print('aaa2:',sl)
#     if not sl:
#         return simplify(0*symbols('x'))
    
#     if not type(sl)=='sympy.core.add.Add':
#         return sl[0]
#     print('aaa:',sl)
    sl1=[parse_expr(item) for item in str(sl[0]).split('+')]
    res=[]
    for k in range(len(sl)-1):
        slt=[parse_expr(item) for item in str(sl[k+1]).split('+')]
        for item1 in sl1:
            for item2 in slt:
                res.append(item1*item2)
        sl1=res
        res=[]
    r=0
    for k in sl1:
        r+=k
    return r

def add_weight(w1,w2):

    r1=get_sum_form(w1)

    r2=get_sum_form(w2)
#     print('ddd:',r1,r2)
    return to_cnf2(r1+r2)

# def add_weight(w1,w2):
#     cofactor=[]
#     cw1=copy.copy(w1)
#     cw2=copy.copy(w2)
# #     print(cw1)
#     for k in cw1:
#         if k in cw2:
#             cofactor.append(k)
#             cw1.remove(k)
#             cw2.remove(k)
#     r1=0
#     r2=0
#     if len(cw1)==1:
#         r1=cw1[0]
#     elif len(cw1)>1:
#         r1=get_sum_form(cw1)
        
#     if len(cw2)==1:
#         r2=cw2[0]
#     elif len(cw2)>1:
#         r2=get_sum_form(cw2)
#     print(w1,w2,cofactor)
#     return cofactor+to_cnf2(r1+r2)

def contract(tdd1,tdd2,key_2_new_key,cont_order,cont_num):
    """The contraction of two TDDs, var_cont is in the form [[4,1],[3,2]]"""
    
    k1=tdd1.node.key
    k2=tdd2.node.key
    w1=tdd1.weight
    w2=tdd2.weight
    
    if k1==-1 and k2==-1:
        if w1==[0]:
            tdd=TDD(tdd1.node)
            tdd.weight=[0]
            return tdd
        if w2==[0]:
            tdd=TDD(tdd1.node)
            tdd.weight=[0]
            return tdd
        tdd=TDD(tdd1.node)
        tdd.weight=mul_weight(w1,w2)
        if cont_num>0:
            tdd.weight[0]*=2**cont_num
#         tdd.weight=simplify(tdd.weight)
        return tdd

    if k1==-1:
        if w1==[0]:
            tdd=TDD(tdd1.node)
            tdd.weight=[0]
            return tdd
        if cont_num ==0 and key_2_new_key[1][k2]==k2:
            tdd=TDD(tdd2.node)
            tdd.weight=mul_weight(w1,w2)
            return tdd
            
    if k2==-1:
        if w2==[0]:
            tdd=TDD(tdd2.node)
            tdd.weight=[0]
            return tdd        
        if cont_num ==0 and key_2_new_key[0][k1]==k1:
            tdd=TDD(tdd1.node)
            tdd.weight=mul_weight(w1,w2)
            return tdd
    
    tdd1.weight=[1]
    tdd2.weight=[1]
    
    temp_key_2_new_key=[]
    temp_key_2_new_key.append(tuple([k for k in key_2_new_key[0][:(k1+1)]]))
    temp_key_2_new_key.append(tuple([k for k in key_2_new_key[1][:(k2+1)]]))
    
    tdd=find_computed_table(['*',tdd1,tdd2,temp_key_2_new_key,cont_num])
    if tdd:
#         tdd.weight=simplify(tdd.weight*w1*w2)
        tdd.weight=mul_weight(tdd.weight,mul_weight(w1,w2))
        tdd1.weight=w1
        tdd2.weight=w2
        return tdd
                
    if cont_order[0][k1]<cont_order[1][k2]:
        the_key=key_2_new_key[0][k1]
        if the_key!='c':
            the_successors=[]
            for k in range(tdd1.node.succ_num):
                res=contract(Slicing(tdd1,k1,k),tdd2,key_2_new_key,cont_order,cont_num)
                the_successors.append(res)
            tdd=normalize(the_key,the_successors)
            insert_2_computed_table(['*',tdd1,tdd2,temp_key_2_new_key,cont_num],tdd)
            tdd.weight=mul_weight(tdd.weight,mul_weight(w1,w2))
        else:
            tdd=TDD(Find_Or_Add_Unique_table(-1))
            tdd.weight=[0]
            for k in range(tdd1.node.succ_num):
                res=contract(Slicing(tdd1,k1,k),tdd2,key_2_new_key,cont_order,cont_num-1)           
                tdd=add(tdd,res)
            insert_2_computed_table(['*',tdd1,tdd2,temp_key_2_new_key,cont_num],tdd)
            tdd.weight=mul_weight(tdd.weight,mul_weight(w1,w2))
    elif cont_order[0][k1]==cont_order[1][k2]:
        the_key=key_2_new_key[0][k1]
        if the_key!='c':
            the_successors=[]
            for k in range(tdd1.node.succ_num):
                res=contract(Slicing(tdd1,k1,k),Slicing(tdd2,k2,k),key_2_new_key,cont_order,cont_num)
                the_successors.append(res)
            tdd=normalize(the_key,the_successors)
            insert_2_computed_table(['*',tdd1,tdd2,temp_key_2_new_key,cont_num],tdd)
            tdd.weight=mul_weight(tdd.weight,mul_weight(w1,w2))
        else:
            tdd=TDD(Find_Or_Add_Unique_table(-1))
            tdd.weight=[0]
            for k in range(tdd1.node.succ_num):
                res=contract(Slicing(tdd1,k1,k),Slicing(tdd2,k2,k),key_2_new_key,cont_order,cont_num-1)           
                tdd=add(tdd,res)
            insert_2_computed_table(['*',tdd1,tdd2,temp_key_2_new_key,cont_num],tdd)
            tdd.weight=mul_weight(tdd.weight,mul_weight(w1,w2))
    else:
        the_key=key_2_new_key[1][k2]
        if the_key!='c':
            the_successors=[]
            for k in range(tdd2.node.succ_num):
                res=contract(tdd1,Slicing(tdd2,k2,k),key_2_new_key,cont_order,cont_num)
                the_successors.append(res)
            tdd=normalize(the_key,the_successors)
            insert_2_computed_table(['*',tdd1,tdd2,temp_key_2_new_key,cont_num],tdd)
            tdd.weight=mul_weight(tdd.weight,mul_weight(w1,w2))
        else:
            tdd=TDD(Find_Or_Add_Unique_table(-1))
            tdd.weight=[0]
            for k in range(tdd2.node.succ_num):
                res=contract(tdd1,Slicing(tdd2,k2,k),key_2_new_key,cont_order,cont_num-1)           
                tdd=add(tdd,res)
            insert_2_computed_table(['*',tdd1,tdd2,temp_key_2_new_key,cont_num],tdd)
            tdd.weight=mul_weight(tdd.weight,mul_weight(w1,w2))
    tdd1.weight=w1
    tdd2.weight=w2
#     tdd.weight=simplify(tdd.weight)
    return tdd
    
def Slicing(tdd,x,c):
    """Slice a TDD with respect to x=c"""

    k=tdd.node.key
    
    if k==-1:
        return tdd.self_copy()
    
    if k<x:
        return tdd.self_copy()
    
    if k==x:
        res=TDD(tdd.node.successor[c])
        res.weight=tdd.node.out_weight[c]
        return res
    else:
        print("Not supported yet!!!")
        
        
def Slicing2(tdd,x,c):
    """Slice a TDD with respect to x=c"""

    k=tdd.node.key
    
    if k==-1:
        return tdd.self_copy()
    
    if k<x:
        return tdd.self_copy()
    
    if k==x:
        res=TDD(tdd.node.successor[c])
#         res.weight=simplify(tdd.node.out_weight[c]*tdd.weight)
        res.weight=mul_weight(tdd.weight,tdd.node.out_weight[c])
        return res
    else:
        print("Not supported yet!!!")        
        
        

def add(tdd1,tdd2):
    """The apply function of two TDDs. Mostly, it is used to do addition here."""
    global global_index_order    
    
    k1=tdd1.node.key
    k2=tdd2.node.key
    
    if tdd1.weight==[0]:
        return tdd2.self_copy()
    
    if tdd2.weight==[0]:
        return tdd1.self_copy()
    
    if tdd1.node==tdd2.node:
        weig=add_weight(tdd1.weight,tdd2.weight)
        if weig==[0]:
            term=Find_Or_Add_Unique_table(-1)
            res=TDD(term)
            res.weight=[0]
            return res
        else:
            res=TDD(tdd1.node)
            res.weight=weig
            return res
        
    if find_computed_table(['+',tdd1,tdd2]):
        return find_computed_table(['+',tdd1,tdd2])
    the_successors=[]
    if k1>k2:
        x=k1
        for k in range(tdd1.node.succ_num):
            res=add(Slicing2(tdd1,x,k),tdd2)
            the_successors.append(res)
    elif k1==k2:
        x=k1
        for k in range(tdd1.node.succ_num):
            res=add(Slicing2(tdd1,x,k),Slicing2(tdd2,x,k))
            the_successors.append(res)        
    else:
        x=k2
        for k in range(tdd2.node.succ_num):
            res=add(tdd1,Slicing2(tdd2,x,k))
            the_successors.append(res)
            
    res = normalize(x,the_successors)
    insert_2_computed_table(['+',tdd1,tdd2],res)
    return res