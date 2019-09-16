import os
import pytest
from rdkit import Chem
from mol_property.similarity import HammingSS, CosineSS
from mol_property.similarity.utils import build_mol_features

cur_dir = os.path.dirname(__file__)

class TestSimilarity():

    smi = "Brc(cc1)cc2c1NCC2"
    mol = Chem.MolFromSmiles(smi)

    zinc_file = os.path.join(cur_dir, "mini_csv.zip")
    feature_file = os.path.join(cur_dir, "fp.arr.npy")
    hamming_index_file = os.path.join(cur_dir, "hamming.index")
    cosine_index_file = os.path.join(cur_dir, "cosine.index")

    @classmethod
    def setup_class(cls):
        build_mol_features(cls.zinc_file, out_file=cls.feature_file)
        assert os.path.exists(cls.feature_file)

    @classmethod
    def teardown_class(cls):
        os.remove(cls.feature_file)
        os.remove(cls.hamming_index_file)
        os.remove(cls.cosine_index_file)

    def test_hammingss(self):
        HammingSS.build_index(feature_file=self.feature_file, 
                              index_file=self.hamming_index_file)
        assert os.path.exists(self.hamming_index_file)
        hs = HammingSS(zinc_file=self.zinc_file, index_file=self.hamming_index_file)
        ret = hs.search_by_mols([Chem.MolFromSmiles(self.smi)])
        assert sorted(ret[0], key=lambda x: x["score"], reverse=True) == ret[0]

    def test_cosiness(self):
        CosineSS.build_index(feature_file=self.feature_file, index_file=self.cosine_index_file)
        cs = CosineSS(zinc_file=self.zinc_file, index_file=self.cosine_index_file)
        ret = cs.search_by_mols([Chem.MolFromSmiles(self.smi)])
        assert sorted(ret[0], key=lambda x: x["score"], reverse=True) == ret[0]

