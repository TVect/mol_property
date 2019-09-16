import os
import numpy as np
import pandas as pd
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs, rdMolDescriptors
import tqdm


def build_mol_features(in_file, out_file):
    df_zinc = pd.read_csv(in_file, compression="zip")
    fp_list = []
    for smi in tqdm.tqdm(df_zinc["smiles"], total=len(df_zinc)):
        tmp_arr = np.array([])
        DataStructs.ConvertToNumpyArray(rdMolDescriptors.GetMACCSKeysFingerprint(Chem.MolFromSmiles(smi)), tmp_arr)
        fp_list.append(tmp_arr)
    fp_arr = np.array(fp_list)
    np.save(out_file, fp_arr)


if __name__ == "__main__":
    build_mol_features("save/standard_csv.zip", "save/fp_arr")
