from qiskit import QuantumCircuit

from .. import DDType
from ..ts import Tensor
from ..tn import TensorNetwork


def cir_2_tn(cir:QuantumCircuit, ddtype:DDType) -> TensorNetwork:
    isPQC = len(cir.parameters) != 0 # PQC: parameterized quantum circuit
    if isPQC:
        if (ddtype is DDType.TDD) or (ddtype is DDType.STDD):
            raise ValueError(f"{ddtype} does not support representing PQC.")

    num_q = cir.num_qubits
    current_indices = []
    
    return TensorNetwork()