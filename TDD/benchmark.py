from qiskit import QuantumCircuit, transpile
from qiskit.circuit import library
import numpy as np


def circuit_generator(qubit_num=1,reps=1,function_name='TwoLocal'
                    ,basis_gates1=['id', 'rz', 'sx', 'x', 'cx']
                    ,basis_gates2=['rx','ry','h','cx']
                    ,error_model = None, error_rate = 0.05
                    ):

    if function_name=='TwoLocal':
        cir = library.TwoLocal(qubit_num, ['ry'],'cx', entanglement='circular', reps=reps)
    if function_name=='ExcitationPreserving':
        cir = library.ExcitationPreserving(qubit_num, mode='fsim', entanglement='full', reps=reps)
    if function_name=='RealAmplitudes':
        cir = library.RealAmplitudes(qubit_num, entanglement='full', reps=reps)
    if function_name=='EfficientSU2':
        cir = library.EfficientSU2(qubit_num, ['rx','h'], entanglement='circular', reps=reps)

    cir1=transpile(cir,basis_gates=basis_gates1)
    cir2=transpile(cir,basis_gates=basis_gates2)

    '''
    add error model start 
    '''
    import random 
    if error_model == 'random_cnot':
        l = cir.num_qubits
        error_rate = error_rate
        error_num = max(int(l*error_rate),1) 
        cnot_list = random.sample(range(l), error_num*2)
        for i in range (error_num):
            cir.cx(cnot_list[2*i],cnot_list[2*i+1])
    if error_model == 'random_flip':
        error_rate = error_rate
        for i in range(len(cir2)):
            if random.random() < error_rate and len(cir2[i].qubits) == 2:
                cir2[i].qubits = (cir2[i].qubits[1], cir2[i].qubits[0])

    if error_model == 'random_shift':
        error_rate = error_rate
        for i in range(len(cir2)):
            gate = cir2[i].operation
            if random.random() < error_rate and gate.is_parameterized():
                gate._params[0] = gate.params[0] + 2 * np.pi * random.random()
    if error_model == 'random_bit_flip':
        def add_error(circuit,error_rate):
            gate_num = len(circuit)
            error_circuit = QuantumCircuit(circuit.num_qubits)
            error_loction = np.random.choice(gate_num,max(int(gate_num*error_rate),1),replace=False)
            count = 0
            for i in circuit:
                error_circuit.append(i)
                if count in error_loction:
                    error_circuit.x(i.qubits[0])
                    # error_circuit.z(i.qubits[0])
                count += 1
            return error_circuit
        cir2 = add_error(cir2,error_rate)
    '''
    add error model end
    '''

    return cir1,cir2,'%s_%i_%i'%(function_name, qubit_num, reps)