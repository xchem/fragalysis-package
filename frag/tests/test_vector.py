from frag.utils.vector_utils import get_exit_vector_for_xe_smi
import unittest


class VectorTest(unittest.TestCase):
    def test_xe_smi(self):
        # Some simple test data - these should all be zero
        input_data = ["[100Xe]c1ccc([101Xe])cc1", "[101Xe]c1ccccc1"]
        output_data = [
            {
                100: 0.0,
                101: 0.0,
            },
            {101: 0.0},
        ]
        # Now iterate over
        for i, smi in enumerate(input_data):
            ans = get_exit_vector_for_xe_smi(smi)
            test_dict = output_data[i]
            for key in test_dict:
                self.assertAlmostEqual(ans[key], test_dict[key], 2)
