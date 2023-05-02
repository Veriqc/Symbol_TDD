from qiskit import QuantumCircuit
from qiskit.circuit import CircuitInstruction
from qiskit.circuit.gate import Gate
from qiskit.circuit.controlledgate import ControlledGate
import qiskit.circuit.library as GateLib
import numpy as np

from collections.abc import Callable
import logging

from .. import DDType, Tensor, Index, TensorNetwork


log = logging.getLogger(__name__)
log.setLevel("DEBUG")

Diagonal_Gates = set([GateLib.HGate, GateLib.IGate, GateLib.PhaseGate,
                      GateLib.RZGate, GateLib.RZZGate, GateLib.SGate, GateLib.SdgGate,
                      GateLib.TGate, GateLib.TdgGate, GateLib.U1Gate, GateLib.ZGate])

Diagonal_Gate_Names = {gate.name for gate in Diagonal_Gates}


def gate_2_data_direct(gate: Gate) -> np.ndarray:
    """
    qiskit gate U_yx reshape to tensor T_yn...y1xn...x1
    ours tensor should be T_x1y1x2y2...xnyn
    """
    try:
        data = gate.to_matrix().reshape((2,) * (2 * gate.num_qubits))
    except:
        # TODO: support qiskit ctrl_state for cases other than full of ones 11...11
        assert isinstance(gate, ControlledGate)
        from scipy.linalg import block_diag
        A = gate.base_gate.to_matrix()
        I = np.eye(*A.shape)
        data = block_diag(*([I] * (2 ** gate.num_ctrl_qubits - 1) + [A]))
        data = data.reshape((2,) * (2 * gate.num_qubits))
        y_p, x_p = np.arange(gate.num_qubits), np.arange(gate.num_qubits, 2 * gate.num_qubits)
        data = np.moveaxis(data, np.hstack((y_p, x_p)), np.hstack((y_p[::-1], x_p[::-1])))

    xpos = np.arange(data.ndim // 2)[::-1] * 2
    ypos = xpos + 1
    data = np.moveaxis(data, np.arange(data.ndim), np.hstack((ypos, xpos)))
    return data

    # is_diagonal = np.count_nonzero(data - np.diag(np.diag(data))) != 0


def instr_2_tensor(gate2data_func: Callable[[Gate], np.ndarray], instr: CircuitInstruction, current_indices: list[Index], cir_qubits: list) -> Tensor:
    """ Create new tensor from the given instruction and update current indices. """
    gate = instr.operation
    data = gate2data_func(gate)

    if isinstance(gate, ControlledGate):
        need_hypridx_list = [True] * gate.num_ctrl_qubits + [gate.base_gate.name in Diagonal_Gate_Names] * (gate.num_qubits - gate.num_ctrl_qubits)
    else:
        need_hypridx_list = [gate.name in Diagonal_Gate_Names] * gate.num_qubits

    qubit_idx_list = [cir_qubits.index(qubit) for qubit in instr.qubits]
    input_indices = [current_indices[q] for q in qubit_idx_list]
    output_indices = [index.create_next(with_hypridx) for index, with_hypridx in zip(input_indices, need_hypridx_list)]
    indices = sum(zip(input_indices, output_indices), ())

    # update current_indices according to output_indices
    for i, q in enumerate(qubit_idx_list):
        current_indices[q] = output_indices[i]

    new_tensor = Tensor(data, indices, gate.name)

    log.debug("instr:%s(%s)", gate.name, ",".join([str(q) for q in qubit_idx_list]))
    log.debug(" => tensor: shape %s, indices %s", data.shape, str(indices))
    # log.debug(" => %s", data)
    return new_tensor


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
        init_tensors = [Tensor(zero_state, (index,)) for index in current_indices]

    # Transform qiskit instructions to tensors
    if ddtype is DDType.STDD:
        raise ValueError("Not implemented")
    elif ddtype is DDType.TrTDD:
        raise ValueError("Not implemented")
    elif ddtype is DDType.ExpTDD:
        raise ValueError("Not implemented")
    else:
        gate2data_func = gate_2_data_direct

    tensors = [instr_2_tensor(gate2data_func, instr, current_indices, cir.qubits) for instr in cir.data]

    # We update lastest index of each qubit to label output 'y' *in-place*
    for q, index in enumerate(current_indices):
        index.update('y', q)
    
    total_tensors = init_tensors + tensors
    log.debug("new TN:")
    for tensor in total_tensors:
        log.debug("%s", str(tensor))
    return TensorNetwork(total_tensors)