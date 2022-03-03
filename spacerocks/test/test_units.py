import unittest
from spacerocks import Units

class TestVector(unittest.TestCase):

    def test_units(self):

        units = Units()
        units.current()

if __name__ == '__main__':
    unittest.main()