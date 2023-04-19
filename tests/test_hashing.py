import unittest

from symtdd.utils import hashing


class TestHashDouble(unittest.TestCase):
    def test_zero(self) -> None:
        self.assertEqual(hashing.hash_double(0.0), 0)
        self.assertEqual(hashing.hash_double(0.5), 0)
        self.assertEqual(hashing.hash_double(0.25), 0)

    def test_005(self) -> None:
        """
        (m=0.8, exp=-4) <= frexp(0.05)
        (2|m|-1)*W = (2*0.8-1) * (2**64)
        """
        self.assertNotEqual(
            hashing.hash_double(0.05), hashing.hash_double(0.05000000000000001)
        )
        self.assertNotEqual(
            hashing.hash_double(0.05), hashing.hash_double(0.04999999999999999)
        )
        self.assertAlmostEqual(hashing.hash_double(0.05), (2 * 0.8 - 1) * (2**64))
