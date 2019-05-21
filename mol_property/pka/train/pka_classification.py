# -*- coding: utf-8
import os
import numpy as np
import xgboost as xgb
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.metrics import classification_report, f1_score, precision_score, recall_score
from sklearn.multiclass import OneVsRestClassifier
from sklearn.externals import joblib

from data_utils import DataUtils

cur_dir = os.path.dirname(__file__)
d_utils = DataUtils(filepath=os.path.join(cur_dir, "data/pKaInWater.csv"))
X_data, y_data_acidic, y_data_basic = d_utils.get_classification_data(feature_type="morgan+macc")
y_data = np.array([y_data_acidic, y_data_basic]).T

# train test split
seed = 7
X_train, X_test, y_train, y_test = train_test_split(X_data, y_data, test_size=0.2, random_state=seed)
print("\n ========================= \n")
print("X_train.shape:", X_train.shape, "X_test.shape", X_test.shape)
print("\n ========================= \n")


def model_evaluation(model, x_input, y_input):
    y_pred = model.predict(x_input)
    print(classification_report(y_true=y_input, y_pred=y_pred))
    print("f1_score [micro]: ", f1_score(y_true=y_input, y_pred=y_pred, average="micro"))
    print("f1_score [macro]: ", f1_score(y_true=y_input, y_pred=y_pred, average="macro"))
    print("precision_score [micro]: ", precision_score(y_true=y_input, y_pred=y_pred, average="micro"))
    print("precision_score [macro]: ", precision_score(y_true=y_input, y_pred=y_pred, average="macro"))
    print("recall_score [micro]: ", recall_score(y_true=y_input, y_pred=y_pred, average="micro"))
    print("recall_score [macro]: ", recall_score(y_true=y_input, y_pred=y_pred, average="macro"))

# sklearn xgboost
xgb_model = xgb.XGBClassifier()
xgb_model_multilabel = OneVsRestClassifier(xgb_model)

# Grid Search CV
# params = {'estimator__max_depth': [4,6,8,10],
#           'estimator__n_estimators': range(80,210,10),
#           'estimator__subsample': [0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
#           'estimator__colsample_bytree': [0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
#           'estimator__eta': [0.01, 0.03, 0.05, 0.08, 0.1, 0.2, 0.3, 0.4, 0.5]}

params = {
          'estimator__colsample_bytree': [0.5],
          'estimator__eta': [0.01],
          'estimator__max_depth': [10],
          # 'estimator__n_estimators': [200],
          'estimator__n_estimators': [80],
          'estimator__subsample': [0.9],
          'estimator__gamma': [0.2]
}

reg = GridSearchCV(xgb_model_multilabel, params, verbose=3, n_jobs=-1, cv=5)
reg.fit(X_train, y_train)
print(reg.best_score_)
print(reg.best_params_)

print("train:", model_evaluation(reg.best_estimator_, X_train, y_train))
print("test:", model_evaluation(reg.best_estimator_, X_test, y_test))

joblib.dump(reg.best_estimator_, "model/pka_classification.pkl")
