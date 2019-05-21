# mol_property

Prediction of pKa from chemical structure using machine learning approaches, and some other property

```
from rdkit import Chem
from mol_property import property_api

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
```
