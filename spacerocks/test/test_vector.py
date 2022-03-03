import unittest
from spacerocks.vector import Vector

class TestVector(unittest.TestCase):

    def test_vector_operations(self):

        a = Vector(1, 1, 1)
        b = Vector(1, 0, 2)

        assert a != b
        assert a + b == Vector(2, 1, 3)
        assert a - b == Vector(0, 1, -1)
        assert 3 * b == Vector(3, 0, 6)
        assert b / 3 == Vector(1/3, 0, 2/3)
        assert a.dot(b) == 1 + 0 + 2
        assert a.cross(b) == Vector(2, -1, -1)
        assert a.norm == 3**0.5
        assert a.unit == Vector(1/3**0.5, 1/3**0.5, 1/3**0.5)

if __name__ == '__main__':
    unittest.main()
