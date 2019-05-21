import pytest
from rdkit import Chem
from mol_property import property_api


class TestPropertyAPI():

    smi = "Brc(cc1)cc2c1NCC2"
    mol = Chem.MolFromSmiles(smi)

    def test_logP(self):
        logP = property_api.get_logP(self.mol)
        assert abs(logP - 2.417) < 0.5
    
    def test_numHBA(self): 
        nHBA = property_api.get_numHBA(self.mol)
        assert nHBA == 1

    def test_nHBD(self):
        nHBD = property_api.get_numHBD(self.mol)
        assert nHBD == 1

    def test_TPSA(self):
        TPSA = property_api.get_polarSurfaceArea(self.mol)
        assert abs(TPSA - 12) < 0.5

    def test_numRotatableBonds(self):
        nRotatableBonds = property_api.get_numRotatableBonds(self.mol)
        assert nRotatableBonds == 0

    def test_pKa(self):
        pKa = property_api.get_pKa(self.mol)
        assert 'basic' in pKa
        assert 'acidic' not in pKa
        assert abs(pKa['basic'] - 3.82) < 0.5

    def test_logS(self):
        logS = property_api.get_logS(self.mol)
        assert abs(logS - (-3.1)) < 0.5

    def test_numRings(self):
        nRings = property_api.get_numRings(self.mol)
        assert nRings == 2

    def test_ruleOfFive(self):
        ruleOfFive = property_api.get_ruleOfFive(self.mol)
        assert not ruleOfFive

    def test_veberRule(self):
        veberRule = property_api.get_veberRule(self.mol)
        assert veberRule

    def test_chemFormula(self):
        chemFormula = property_api.get_chemicalFormula(self.mol)
        assert chemFormula == "C8H8BrN"

    def test_Mass(self):
        mass = property_api.get_molecularMass(self.mol)
        assert abs(mass - 196.98) < 5

    def test_SMILES(self):
        smiles = property_api.get_SMILES(self.mol)
        
