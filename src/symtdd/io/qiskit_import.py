from qiskit import QuantumCircuit
from qiskit.circuit import CircuitInstruction
from qiskit.circuit.gate import Gate
import numpy as np

from collections.abc import Callable

from .. import DDType, Tensor, Index, TensorNetwork


def gate_2_data_direct(gate: Gate) -> np.ndarray:
    """
    qiskit gate U_yx reshape to tensor T_yn...y1xn...x1
    ours tensor should be T_x1y1x2y2...xnyn
    """
    data = gate.to_matrix().reshape((2,) * gate.num_qubits)
    xpos = np.arange(data.ndim // 2)[::-1] * 2
    ypos = xpos + 1
    data = np.moveaxis(data, np.arange(data.ndim), np.hstack((ypos, xpos)))
    return data


def instr_2_tensor(gate2data_func: Callable[[Gate], np.ndarray], instr: CircuitInstruction, input_indices: list) -> Tensor:
    pass


def cir_2_tn(cir:QuantumCircuit, ddtype:DDType) -> TensorNetwork:
    if isPQC := (len(cir.parameters) != 0): # PQC: parameterized quantum circuit
        if (ddtype is DDType.TDD) or (ddtype is DDType.STDD):
            raise ValueError(f"{ddtype} does not support representing PQC.")

    num_q = cir.num_qubits
    current_indices = [Index('x', q) for q in range(num_q)]

    # Construct initial state tensors
    if ddtype is DDType.STDD:
        raise ValueError("Not implemented")
    else:
        zero_state = np.array([1, 0])
        init_tensors = [Tensor(zero_state, index) for index in current_indices]

    # Transform qiskit instructions to tensors
    if ddtype is DDType.STDD:
        raise ValueError("Not implemented")
    elif ddtype is DDType.TrTDD:
        raise ValueError("Not implemented")
    elif ddtype is DDType.ExpTDD:
        raise ValueError("Not implemented")
    else:
        gate2data_func = gate_2_data_direct

    tensors = [instr_2_tensor(gate2data_func, instr, current_indices) for instr in cir.data]

    # We update lastest index of each qubit to label output 'y' *in-place*
    for q, index in enumerate(current_indices):
        index.update('y', q)
    
    return TensorNetwork(tensors)