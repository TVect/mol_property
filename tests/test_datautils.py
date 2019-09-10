import pytest
from mol_property.pka import data_utils

class TestDatautils():
    
    d_utils = data_utils.DataUtils()

    def test_get_classification_data(self):
        ret = self.d_utils.get_classification_data(feature_type="morgan")
        assert len(ret) == 3
        assert ret[0].shape == (7466, 1024)
        assert (ret[0].shape, ret[1].shape, ret[2].shape) == ((7466, 1024), (7466,), (7466,))
        
        ret = self.d_utils.get_classification_data(feature_type="macc")
        assert len(ret)  == 3 
        assert (ret[0].shape, ret[1].shape, ret[2].shape) == ((7466, 167), (7466,), (7466,))

        ret = self.d_utils.get_classification_data(feature_type="morgan+macc")
        assert len(ret)  == 3 
        assert (ret[0].shape, ret[1].shape, ret[2].shape) == ((7466, 1191), (7466,), (7466,))

    def test_get_regression_data(self):
        ret = self.d_utils.get_regression_data(data_category="all")
        assert len(ret) == 2
        assert (len(ret[0]), len(ret[1])) == (7911, 7911)

        ret = self.d_utils.get_regression_data(data_category="acidic_only")
        assert len(ret) == 2
        assert (len(ret[0]), len(ret[1])) == (3614, 3614)

        ret = self.d_utils.get_regression_data(data_category="basic_only")
        assert len(ret) == 2
        assert (len(ret[0]), len(ret[1])) == (4297, 4297)

