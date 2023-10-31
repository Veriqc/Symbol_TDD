import numpy as np
from TDD.TDD import Index,Ini_TDD
from TDD.TDD_Q import cir_2_tn
from TDD.TN import Tensor

from qiskit import QuantumCircuit
from sympy import *


import time, csv

def simulation(input_cir, Benchmark_Name=None,symbolic=True, unique_table_reset=True,output_file=None):
    if isinstance(input_cir,QuantumCircuit):
        cir=input_cir
    else:
        cir=QuantumCircuit.from_qasm_file(input_cir)
        if not Benchmark_Name:
            Benchmark_Name=input_cir
    tn,indices=cir_2_tn(cir)
    if symbolic:
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
    tdd,Max_node_num=tn.cont(max_node=True)
    Time=time.time()-t_start
    from TDD.TDD import get_unique_table_num as SymTDD_gu1
    from TDD.SymTDD.BDD import get_unique_table_num as SymTDD_gu2
    output_dict = {'Benchmark Name':Benchmark_Name,
            'Qubit num.':tn.qubits_num,
            'Gate num.':len(cir.data),
            'Time':Time,
            'Node num. max':Max_node_num,
            'Node num. final':tdd.node_number(),
            'gu1': SymTDD_gu1(),
            'gu2': SymTDD_gu2()}

    return tdd, output_dict

'''
Inverse test
'''
def Inverse_test(input_file,optimizer=None,output_file=None):
    cir=QuantumCircuit.from_qasm_file(input_file)
    cir2=cir.inverse()
    cir3=cir.compose(cir2)
    tn,indices=cir_2_tn(cir3)
    for k in range(tn.qubits_num):
        x_k='x'+str(k)
        xn_k='xn'+str(k)
        s=Symbol(x_k)
        ns=Symbol(xn_k)
        U=np.array([ns,s])
        # U=np.array([0,1])
        temp_ts=Tensor(U,[Index(x_k)])
        tn.tensors.insert(0,temp_ts)
        tn.tensors.append(Tensor(np.array([[ns],[s]]),[Index('y%i'%k)]))
        if not x_k in indices:
            indices.append(x_k)
    t_start=time.time()
    Ini_TDD(indices,n=300,type='SymTDD')
    tdd,Max_node_num=tn.cont(optimizer=optimizer,max_node=True)

    Time=time.time()-t_start
    # print('Benchmark Name:',input_file)
    # print('Time:',Time)
    # print('Qubit num.:',tn.qubits_num)
    # print('Gate num.:',len(cir3.data))
    # print('Node num. max:',Max_nodes[0])
    # print('Node num. final:',tdd.node_number())
    # print('gu1:',gu1())
    # print('gu2:',gu2())
    if output_file:
        with open(output_file, 'a', newline='') as csvfile:
            fieldnames = ['Benchmark Name', 'Qubit num.','Gate num.','Time','Node num. max','Node num. final','gu1','gu2']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writerow({'Benchmark Name':input_file,
            'Qubit num.':tn.qubits_num,
            'Gate num.':len(cir3.data),
            'Time':Time,
            'Node num. max':Max_node_num,
            'Node num. final':tdd.node_number(),
            'gu1':gu1(),
            'gu2':gu2()})
    return tdd