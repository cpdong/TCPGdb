#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 20 10:39:52 2024

@author: cpdong
"""

import argparse;
import pandas as pd
import numpy as np
import pickle, time
from pathlib import Path
from sklearn.metrics import roc_auc_score
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

def Grid_Search_CV_LASSO(x_train, y_train, cv_fold=5):
    from sklearn.model_selection import GridSearchCV
    from sklearn.linear_model import Lasso
    from sklearn.model_selection import StratifiedKFold
    from sklearn.model_selection import LeaveOneOut

    param_grid = [{"max_iter": [100],
                   "alpha": np.linspace(start=0, stop=1, num=101) }]
    
    estimator = Lasso()
    grid_search = GridSearchCV(estimator = estimator,
                               param_grid = param_grid,
                               cv = LeaveOneOut(),
                               scoring='neg_mean_absolute_error');
    
    grid_search.fit(x_train,y_train)
    best_params = grid_search.best_params_
    grid_search_cv_results = pd.DataFrame(grid_search.cv_results_)
    best_estimator = Lasso().set_params(**best_params)

    return best_params, best_estimator, grid_search_cv_results


if __name__ == '__main__':
    cpus_num = 4
    cv_fold = 5
    os.chdir('/path/to/workdir')

    geneList = [l.rstrip() for l in open('./CD8CD4_top50list.txt')]
    features = list(set(geneList))
    
    trainData = pd.read_csv('./GSE223655_data.tsv', header=0, index_col=0, sep='\t').T
    y_train = pd.factorize(trainData['group'], sort=True)[0]
    x_train = trainData[features]
    
    testData = pd.read_csv('./GSE160154_data.tsv', header=0, index_col=0, sep='\t').T

    best_params, best_estimator, grid_search_cv_results = Grid_Search_CV_LASSO(x_train, y_train, cv_fold=cv_fold)
    best_estimator.fit(x_train, y_train);

    test_pred = best_estimator.predict(testData[features])
    test_auc = roc_auc_score(testData['group'], test_pred)
    
    out_df = pd.DataFrame({'pred': test_pred, 'label':testData['group'] })
    out_df.to_csv('./GSE160154_LASSO_pred.tsv', sep='\t')
#        