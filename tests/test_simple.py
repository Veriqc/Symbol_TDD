import unittest
from qiskit.circuit.random import random_circuit
# from qiskit import Aer
from qiskit_aer import AerSimulator
from qiskit import transpile
import numpy as np

from symtdd.io.qiskit_import import cir_2_tn
from symtdd import DDType


def test_null() -> None:
    pass


class TestTN(unittest.TestCase):
    def setUp(self) -> None:
        self.cir = random_circuit(5, 100, measure=False, seed=12345)
        cir = self.cir.copy()
        cir.save_statevector()
        # backend = Aer.get_backend('aer_simulator_statevector')
        backend = AerSimulator(method='statevector')
        cir2 = transpile(cir, backend)
        job = backend.run(cir2)
        result = job.result()
        self.gold_data = result.get_statevector()
        return super().setUp()
    
    def test_ts(self) -> None:
        tn = cir_2_tn(self.cir, DDType.TS)
        ts = tn.contract()
        ts.sort()
        data = ts.data
        orders = np.arange(data.ndim)
        data_qiskit_order = np.moveaxis(data, orders, orders[::-1])
        self.assertTrue(np.allclose(self.gold_data, data_qiskit_order.flatten()))