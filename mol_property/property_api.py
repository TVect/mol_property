# -*- coding: utf-8 -*-

import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs

from rdkit.Chem.Descriptors import ExactMolWt    # Molecular Weight
from rdkit.Chem.Crippen import MolLogP    # xlogP
# Net Charge
from rdkit.Chem.rdMolDescriptors import CalcNumRotatableBonds    # Rotate Bonds
from rdkit.Chem.rdMolDescriptors import CalcTPSA    # (topology?) Polar Surface Area
from rdkit.Chem.rdMolDescriptors import CalcNumHBD    # Hydrogen Donors
from rdkit.Chem.rdMolDescriptors import CalcNumHBA    # Hydrogen Acceptors
from rdkit.Chem.rdMolDescriptors import CalcNumRings    # number of rings for a molecule
from rdkit.Chem.rdMolDescriptors import CalcMolFormula    # mol formula
from .solubility.esol import ESOLCalculator
from .pka.predictor import PkaPredictor

esol_calculator = ESOLCalculator()
pka_predictor = PkaPredictor()

def get_logP(mol):
    ''' clogP 或 LogP '''
    return MolLogP(mol)

def get_logD(mol):
    ''' logD '''
    raise NotImplmentedError

def get_physioCharge(mol):
    ''' Net Charge 或Formal Charge 或 Physiological Charge '''
    raise NotImplmentedError

def get_numHBA(mol):
    ''' Hydrogen Acceptor Count '''
    return CalcNumHBA(mol)

def get_numHBD(mol):
    ''' CalcNumHBD '''
    return CalcNumHBD(mol)

def get_polarSurfaceArea(mol):
    ''' Polar Surface Area '''
    return CalcTPSA(mol)

def get_numRotatableBonds(mol):
    ''' Rotatable Bond Count '''
    return CalcNumRotatableBonds(mol)

def get_pKa(mol):
    '''pKa (Strongest Acidic) or pKa (Strongest Basic)'''
    return pka_predictor.predict([mol])[0]

def get_logS(mol):
    ''' Water Solubility 或 logS '''
    return esol_calculator.calc_esol(mol)

def get_numRings(mol):
    ''' Number of Rings '''
    return CalcNumRings(mol)

def get_bioavailability(mol):
    ''' Bioavailability '''
    raise NotImplmentedError

def get_ruleOfFive(mol):
    ''' Rule of Five 
    相对分子质量小于等于500；脂水分配系数小于等于5；氢键供体小于等于5；氢键受体小于10
    '''
    if (get_molecularMass(mol) > 500) and (get_logP(mol) <=5) and (get_numHBD(mol) <= 5) and (get_numHBA(mol) < 10):
        return True
    return False

def get_veberRule(mol):
    ''' Veber's Rule 
    可旋转键数小于等于10或者分子中环数目小于等于4；极性表面积或者氢键总数小于等于12
    '''
    if ((get_numRotatableBonds(mol) <= 10) or (get_numRings(mol) <=4)) and \
       ((get_polarSurfaceArea(mol) <= 12) or (get_numHBD(mol) + get_numHBA(mol) <= 12)):
        return True
    return False

def get_chemicalFormula(mol):
    ''' Chemical Formula '''
    return CalcMolFormula(mol)

def get_molecularMass(mol):
    ''' Molecular Mass '''
    return ExactMolWt(mol)

def get_SMILES(mol):
    ''' SMILES '''
    return Chem.MolToSmiles(mol)


if __name__ == "__main__":
    smi = "Brc(cc1)cc2c1NCC2"
    mol = Chem.MolFromSmiles(smi)

    print("logP:", get_logP(mol))
    # print("logD:", get_logD(mol))
    # print("physioCharge:", get_physioCharge(mol))
    print("HBA:", get_numHBA(mol))
    print("HBD:", get_numHBD(mol))
    print("TPSA:", get_polarSurfaceArea(mol))
    print("Rotatable Bond:", get_numRotatableBonds(mol))
    print("pKa:", get_pKa(mol))
    print("logS:", get_logS(mol))
    print("numRings:", get_numRings(mol))
    # print("Bioavailability", get_bioavailability(mol))
    print("ruleOfFive:", get_ruleOfFive(mol))
    print("veberRule:", get_veberRule(mol))
    print("ChemicalFormula:", get_chemicalFormula(mol))
    print("Mass:", get_molecularMass(mol))
    print("smiles:", get_SMILES(mol))
