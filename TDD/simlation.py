import time
import numpy as np
from qiskit import QuantumCircuit

from TDD.TDD import Index,Ini_TDD
from TDD.TDD_Q import cir_2_tn
from TDD.TN import Tensor

def SymTDD_simulation(input_cir,symbolic=True, unique_table_reset=True, output_file=None):
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
    tdd,Max_node_num=tn.cont(max_node=True)
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
