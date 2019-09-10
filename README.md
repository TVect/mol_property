# mol_property

[![Build Status](https://api.travis-ci.com/TVect/mol_property.svg?branch=master)](https://api.travis-ci.com/TVect/mol_property.svg)    [![codecov](https://codecov.io/gh/TVect/mol_property/branch/master/graph/badge.svg)](https://codecov.io/gh/TVect/mol_property)

Prediction of pKa from chemical structure using machine learning approaches.

The repository also calculate some other property, such as solubility, number of Rotatable Bonds.

## Properties

- **pKa**

ref: [Prediction of pKa from chemical structure using free and open-source tools](https://cfpub.epa.gov/si/si_public_file_download.cfm?p_download_id=535243&Lab=NCCT)

- **logS**

ref: [github: solubility](https://github.com/PatWalters/solubility.git)

## install

- **python package**

```
git clone https://github.com/TVect/mol_property.git
cd mol_property
python setup.py install
```

- **conda package**

reference: https://github.com/fastai/fastai/blob/master/conda/meta.yaml

```
git clone https://github.com/TVect/mol_property.git
cd mol_property
python setup.py sdist
conda build conda/
```


## Usage

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

