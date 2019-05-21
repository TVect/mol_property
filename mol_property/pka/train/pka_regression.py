# -*- coding: utf-8
import os
import numpy as np
import xgboost as xgb
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.metrics import r2_score, mean_squared_error
from sklearn.externals import joblib

from data_utils import DataUtils

# basic_only | acidic only
DATA_CATEGORY = "basic_only"
FEATURE_TYPE = "morgan+macc"

cur_dir = os.path.dirname(__file__)
d_utils = DataUtils(filepath=os.path.join(cur_dir, "data/pKaInWater.csv"))
X_data, y_data = d_utils.get_regression_data(data_category=DATA_CATEGORY, feature_type=FEATURE_TYPE)

# train test split
seed = 7
X_train, X_test, y_train, y_test = train_test_split(X_data, y_data, test_size=0.2, random_state=seed)
print("\n ========================= \n")
print("X_train.shape:", X_train.shape, "X_test.shape", X_test.shape)
print("\n ========================= \n")


def model_evaluation(model, x_input, y_input):
    y_pred = model.predict(x_input)
    rmse_value = np.sqrt(mean_squared_error(y_true=y_input, y_pred=y_pred))
    r2_value = r2_score(y_true=y_input, y_pred=y_pred)
    return rmse_value, r2_value


# sklearn xgboost
xgb_model = xgb.XGBRegressor(tree_method="hist")

# Grid Search CV
params = {'max_depth': [4,6,8,10],
          'n_estimators': range(80,210,10),
          'subsample': [0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
          'colsample_bytree': [0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
          'eta' : [0.01, 0.03, 0.05, 0.08, 0.1, 0.2, 0.3, 0.4, 0.5]}

params = {'colsample_bytree': [0.5], 'eta': [0.01], 'max_depth': [10], 'n_estimators': [420], 'subsample': [0.9],
          'gamma': [0.2]}

reg = GridSearchCV(xgb_model, params, verbose=3, n_jobs=-1, cv=5)
reg.fit(X_train, y_train)
print(reg.best_score_)
print(reg.best_params_)

print("train:", model_evaluation(reg.best_estimator_, X_train, y_train))
print("test:", model_evaluation(reg.best_estimator_, X_test, y_test))

joblib.dump(reg.best_estimator_, "model/pka_{}_regression.pkl".format(DATA_CATEGORY.split("_")[0]))
