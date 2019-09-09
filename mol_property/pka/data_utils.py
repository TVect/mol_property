# -*- coding: utf-8 -*-

import os
import numpy as np
import pandas as pd
import rdkit
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Descriptors

DATA_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "./train/data/pKaInWater.csv")

def rdkit_numpy_convert(fp):
    output = []
    for f in fp:
        arr = np.zeros((1,))
        DataStructs.ConvertToNumpyArray(f, arr)
        output.append(arr)
    return np.asarray(output)


class DataUtils(object):

    def __init__(self, filepath=DATA_PATH):
        self.filepath = filepath

        self.df_pka = pd.read_csv(self.filepath)
        self.df_pka_acidic = self.df_pka[self.df_pka["basicOrAcidic"] == "acidic"]
        self.df_pka_basic = self.df_pka[self.df_pka["basicOrAcidic"] == "basic"]

    def describe(self):
        print("Unique: {} / {}".format(len(self.df_pka["Smiles"].unique()),
                                       self.df_pka.shape[0]))
        print("basic Unique: {} / {}".format(len(self.df_pka_basic["Smiles"].unique()),
                                             self.df_pka_basic.shape[0]))
        acidic_only_cnt = len(set(self.df_pka_acidic["Smiles"].unique()) - set(self.df_pka_basic["Smiles"].unique()))
        basic_only_cnt = len(set(self.df_pka_basic["Smiles"].unique()) - set(self.df_pka_acidic["Smiles"].unique()))
        both_pka_cnt = len(set(self.df_pka_basic["Smiles"].unique()) & set(self.df_pka_acidic["Smiles"].unique()))
        print("acidic_only_cnt: {}, basic_only_cnt: {}, both_pka_cnt: {}".format(acidic_only_cnt, basic_only_cnt,
                                                                                 both_pka_cnt))

    def get_regression_data(self, data_category="all", feature_type="morgan"):
        '''
        :param data_category: all | acidic_only | basic_only
        :param type: morgan | macc | morgan+macc
        :return:
        '''
        df_tmp = self.df_pka
        if data_category == "basic_only":
            df_tmp = self.df_pka_basic
        elif data_category == "acidic_only":
            df_tmp = self.df_pka_acidic

        mols = []
        targets = []
        for row in df_tmp[["pKa", "Smiles"]].iterrows():
            pka, smi = row[1]
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                print(smi)
            else:
                mols.append(mol)
                targets.append(pka)
        return self.get_molecular_features(mols, feature_type), targets

    def get_classification_data(self, feature_type="morgan+macc"):
        '''
        :param type: morgan | macc | morgan+macc
        :return:
        '''
        smi_dict = {}
        for row in self.df_pka[["basicOrAcidic", "Smiles"]].iterrows():
            basicOrAcidic, smi = row[1]
            if Chem.MolFromSmiles(smi) is not None:
                if smi not in smi_dict:
                    smi_dict[smi] = {"basic": 0, "acidic": 0}
                smi_dict[smi][basicOrAcidic] = 1
        df_smi = pd.DataFrame(smi_dict).transpose()
        return self.get_molecular_features([Chem.MolFromSmiles(mol) for mol in df_smi.index], feature_type), \
               df_smi["acidic"].values, df_smi["basic"].values

    @staticmethod
    def get_molecular_features(mols, feature_type="morgan+macc"):
        '''
        :param mols: moleculars
        :param feature_type:
        :return: morgan | macc | morgan+macc
        '''
        fp_morgan = [AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024) for mol in mols]
        fp_macc = [AllChem.GetMACCSKeysFingerprint(mol) for mol in mols]
        if feature_type == "morgan":
            return rdkit_numpy_convert(fp_morgan)
        elif feature_type == "macc":
            return rdkit_numpy_convert(fp_macc)
        elif feature_type == "morgan+macc":
            return np.concatenate([rdkit_numpy_convert(fp_morgan), rdkit_numpy_convert(fp_macc)], axis=1)

