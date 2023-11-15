import time
import numpy as np
from qiskit import QuantumCircuit

from TDD.TDD import Index,Ini_TDD
from TDD.TDD_Q import cir_2_tn, get_real_qubit_num
from TDD.TN import Tensor

def SymTDD_simulation(input_cir,symbolic=True, optimizer=None, unique_table_reset=True, output_file=None):
    if isinstance(input_cir,QuantumCircuit):
        cir=input_cir
    else:
        cir=QuantumCircuit.from_qasm_file(input_cir)
    tn,indices=cir_2_tn(cir)
    if symbolic:
        from sympy import Symbol
        for k in range(tn.qubits_num):
            x_k='x'+str(k)
            xn_k='xn'+str(k)
            s=Symbol(x_k)
            ns=Symbol(xn_k)
            U=np.array([ns,s])
            # U=np.array([0,1])
            temp_ts=Tensor(U,[Index(x_k)])
            tn.tensors.insert(0,temp_ts)
            if not x_k in indices:
                indices.append(x_k)
    t_start=time.time()
    Ini_TDD(indices,n=300,type='SymTDD',unique_table_reset=unique_table_reset)
    # Ini_TDD(indices,n=300,unique_table_reset=unique_table_reset)
    tdd,Max_node_num=tn.cont(optimizer = optimizer, max_node=True)
    Time=time.time()-t_start
    from TDD.TDD import get_unique_table_num as SymTDD_gu1
    from TDD.SymTDD.BDD import get_unique_table_num as SymTDD_gu2
    output_dict = {
            'Qubit num.':tn.qubits_num,
            'Gate num.':len(cir.data),
            'Time':Time,
            'Node num. max':Max_node_num,
            'Node num. final':tdd.node_number(),
            'gu1': SymTDD_gu1(),
            'gu2': SymTDD_gu2()}

    return tdd, output_dict

def TrDD_simulation(cir, Benchmark_Name=None, unique_table_reset=True, add_inputs_list=None, optimizer=None):

    n=get_real_qubit_num(cir)
    tn,indices=cir_2_tn(cir)
    Parameter_num=len(cir.parameters)
    #add y indices
    indices2=[]
    for i, item in enumerate(indices):
        indices2.append(item)
        if item[0]=='y':
            num=int(item.replace('y',''))
            # indices2.append('z%i'%num)
    t_start=time.time()
    #add sin cos indices
    '''
    TrDD process start
    '''
    sin_str=set()
    cos_str=set()
    for tensor in tn.tensors:
        for element in tensor.data.flatten(): 
            from sympy.core.expr import Expr 
            if isinstance(element,Expr):
                for symbol in element.free_symbols:
                    sin_str.add('sin('+str(symbol)+')')
                    cos_str.add('cos('+str(symbol)+')')
    sin_str=list(sin_str)
    sin_str.sort()
    cos_str=list(cos_str)
    cos_str.sort()
    sym_str=[]
    for i in range(len(sin_str)):
        sym_str.append(sin_str[i])
        sym_str.append(cos_str[i])

    # TDD process
    Ini_TDD(indices2,sym_str,type='TrDD',unique_table_reset=unique_table_reset)

    '''
    TrDD process end
    '''

    if add_inputs_list:
        from TDD.TDD_Q import add_inputs
        add_inputs(tn,add_inputs_list)
    #How long start cont
    start_cont=time.time()-t_start
    tdd, Max_node_num=tn.cont(optimizer=optimizer,max_node=True)
    #How long cont process
    cont_time=time.time()-start_cont-t_start

    from TDD.TDD import get_unique_table_num as gu1
    from TDD.TrDD.BDD import get_unique_table_num as gu2
    output_dict={'Benchmark Name':Benchmark_Name,
            'Parameter num': Parameter_num,
            'Qubit num.':tn.qubits_num,
            'Gate num.':len(cir.data),
            'Time':cont_time,
            'Node num. max':Max_node_num,
            'Node num. final':tdd.node_number(),
            'gu1':gu1(),
            'gu2':gu2()
    }
    return tdd, output_dict

def TrDD_verify(cir, cir2, Benchmark_Name=None, unique_table_reset=True, add_inputs_list=None, optimizer=None):

    n=get_real_qubit_num(cir)
    cir_composed  = cir.compose(cir2.inverse())

    tensor_count=0
    Parameter_tensor_location={}

    from sympy import IndexedBase, Symbol
    for gate in cir_composed.data:
        if gate.operation.is_parameterized():
            param_expr = gate.operation.params[0] 

            sym_str=param_expr.parameters.copy().pop().name
            def get_numbers(s):
                import re
                # 使用正則表達式找到所有匹配"\d+"的子串，即連續的一個或多個數字
                numbers = re.findall("\d+", s)
                return int(numbers[0])
            if '[' in sym_str: #利用IndexedBase去對應qiskit生成的Parameter(θ[0])的變數
                sym_str=sym_str.replace("[","").replace("]","")
                sym_order=get_numbers(sym_str)
                sym_str=sym_str.replace(str(sym_order),"")
                sym_base = IndexedBase(sym_str)
                s=sym_base[sym_order]
            else:
                s=Symbol(sym_str)
            Parameter_tensor_location[tensor_count]= s
        tensor_count += 1
    tn, indices = cir_2_tn(cir_composed)
    Parameter_num = cir.num_parameters
    t_start=time.time()

    '''
    TrDD process start
    '''
    # add sin cos indices
    sin_str=set()
    cos_str=set()
    
    for tensor in tn.tensors:
        for element in tensor.data.flatten(): 
            from sympy.core.expr import Expr 
            if isinstance(element,Expr):
                for symbol in element.free_symbols:
                    sin_str.add('sin('+str(symbol)+')')
                    cos_str.add('cos('+str(symbol)+')')

    sin_str=list(sin_str)
    sin_str.sort()
    cos_str=list(cos_str)
    cos_str.sort()
    sym_str=[]
    for i in range(len(sin_str)):
        sym_str.append(sin_str[i])
        sym_str.append(cos_str[i])

    # TDD process
    Ini_TDD(indices,sym_str,type='TrDD',unique_table_reset=unique_table_reset)

    '''
    TrDD process end
    '''

    if add_inputs_list:
        from TDD.TDD_Q import add_inputs
        add_inputs(tn,add_inputs_list)

    def cont2(tn,l1):
        from TDD.TDD import cont, get_identity_tdd
        max_node_num=0
        tdd=get_identity_tdd()
        
        tensor_location1= l1 -1
        tensor_location2= l1
        tensors_len = len(tn.tensors)
        Parameter_location = Parameter_tensor_location.keys()

        #Cont between two parameters(contain).
        def inner_cont (tdd, tensor_location1, tensor_location2, max_node_num):
            have_cont = False
            if tensor_location1 not in Parameter_location and tensor_location1>=0:
                tdd1=tn.tensors[tensor_location1].tdd()
                tdd=cont(tdd1,tdd)
                max_node_num=max(max_node_num,tdd.node_number())
                tensor_location1 -= 1
                have_cont = True
            if tensor_location2 not in Parameter_location and tensor_location2 < tensors_len:
                tdd2=tn.tensors[tensor_location2].tdd()
                tdd=cont(tdd,tdd2)
                max_node_num=max(max_node_num,tdd.node_number())
                tensor_location2 += 1
                have_cont = True
            if have_cont:
                tdd, tensor_location1, tensor_location2, max_node_num = inner_cont(tdd, tensor_location1, tensor_location2, max_node_num)
            else:
                if tensor_location1>=0:
                    tdd1=tn.tensors[tensor_location1].tdd()
                    tdd=cont(tdd1,tdd)
                    max_node_num=max(max_node_num,tdd.node_number())
                if tensor_location2 < tensors_len:
                    tdd2=tn.tensors[tensor_location2].tdd()
                    tdd=cont(tdd,tdd2)
                    max_node_num=max(max_node_num,tdd.node_number())

                tensor_location1 -=1
                tensor_location2 +=1
            return tdd, tensor_location1, tensor_location2, max_node_num
        
        #check whether parameters exist in the tdd
        from functools import lru_cache
        @lru_cache(maxsize=None)
        def parameter_inorder_check(node):
            #inorder triverse the tdd
            if node.key == -1:
                return True
            #left node
            if node.out_weight[0].node.key !=-1:
                return False
            #right node
            if node.out_weight[1].node.key !=-1:
                return False

            #check left node
            if not parameter_inorder_check(node.successor[0]):
                return False
            #check right node
            if not parameter_inorder_check(node.successor[1]):
                return False
            
            return True

        while tensor_location1 >= 0 or tensor_location2 < len(tn.tensors):
            # print('tensor_location',tensor_location1,tensor_location2)
            
            tdd, tensor_location1, tensor_location2, max_node_num = inner_cont(tdd, tensor_location1, tensor_location2, max_node_num)        
            if not parameter_inorder_check(tdd.node):
                break
            
        return tdd, max_node_num


    l1=len(cir.data)
    #How long start cont
    start_cont=time.time()-t_start
    tdd, Max_node_num= cont2(tn,l1)
    #How long cont process
    cont_time=time.time()-start_cont-t_start
    
    identity_cir = QuantumCircuit(n)
    for i in range(n):
        identity_cir.id(i)
    identity, ouput_dict = TrDD_simulation(identity_cir, unique_table_reset= False)
    if tdd.node !=  identity.node:
        equivalent = 'not_equivalent'
    elif tdd.weight == identity.weight:
        equivalent = 'equivalent'
    else:
        equivalent = 'equivalent_up_to_global_phase'

    from TDD.TDD import get_unique_table_num as gu1
    from TDD.TrDD.BDD import get_unique_table_num as gu2
    output_dict={'Benchmark Name':Benchmark_Name,
            'Parameter num.': Parameter_num,
            'Qubit num.':tn.qubits_num,
            'Gate num1.':len(cir.data),
            'Gate num2.':len(cir2.data),
            'Node num. max':Max_node_num,
            'Node num. final':tdd.node_number(),
            'gu1':gu1(),
            'gu2':gu2(),
            'Equivalent':equivalent , 
            'Time':cont_time

    }
    return tdd, output_dict