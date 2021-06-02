#!/usr/bin/env python
# coding: utf-8

# In[111]:


# New approach lets attack building the shell.
import glob
import io
import math
import os
import re
import time
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tensorflow as tf
from lightgbm import LGBMClassifier
from matplotlib import pyplot
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.feature_selection import SelectFromModel, SelectKBest, chi2
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import explained_variance_score
from sklearn.model_selection import cross_val_score, train_test_split
from sklearn.preprocessing import MinMaxScaler, OneHotEncoder
from tensorflow import keras
from tensorflow.keras import Sequential, backend, layers
from tensorflow.keras.layers import Dense, Dropout


# In[103]:


# make classes that expect an input of certain input For DNA, RNA, and methylation


class molecule_preprocessing:
    def __init__(self, path):
        self.path = path

        pass

    # pass as a parameter following testing
    migrated_cgIDs_list = pd.read_csv(
        "/home/jovyan/storage/Machine_Learning/infinium-methylationepic-vs-humanmethylation450.csv"
    )

    def collect_migrated_cgIDs(df):
        """
        read the file into a numpy array that can be passed to the reduce function
        """
        migrated_cgIDs = df.values.flatten()

        return migrated_cgIDs

    def reduce(df, migrated_cgIDs):
        dat = df[np.intersect1d(df.columns, migrated_cgIDs)]

        return dat

    def read_cut(self, file, start, end, molecule):

        if molecule == "RNA":
            sampleIDs = pd.read_table(file, nrows=0, delimiter=",").columns
            df = (
                pd.read_table(
                    file,
                    delimiter=",",
                    skiprows=start,
                    nrows=end - start,
                    names=sampleIDs,
                )
                .fillna(0)
                .set_index("Ensembl")
                .T
            )
            if "__no_feature" in df.columns:
                df = df.drop("__no_feature", axis=1)
            if "__ambiguous" in df.columns:
                df = df.drop("__ambiguous", axis=1)
            if "__too_Low_aQual" in df.columns:
                df = df.drop("__too_low_aQual", axis=1)
            if "__not_aligned" in df.columns:
                df = df.drop("__not_aligned", axis=1)
            if "__alignment_not_unique" in df.columns:
                df = df.drop("__alignment_not_unique", axis=1)

            df.rename(columns={"Ensembl": "samples"}, inplace=True)

        if molecule == "DNA":
            sampleIDs = pd.read_table(
                file, delimiter="\t", compression="gzip", nrows=0
            ).columns
            df = (
                pd.read_table(
                    file,
                    delimiter="\t",
                    compression="gzip",
                    skiprows=start,
                    nrows=end - start,
                    names=sampleIDs,
                )
                .fillna(0)
                .set_index("sample")
                .T
            )

        if molecule == "methylation":
            sampleIDs = pd.read_table(
                file, delimiter="\t", compression="gzip", nrows=0
            ).columns
            df = (
                pd.read_table(
                    file,
                    delimiter="\t",
                    compression="gzip",
                    skiprows=start,
                    nrows=end - start,
                    names=sampleIDs,
                )
                .fillna(0)
                .set_index("Composite Element REF")
                .T
            )
            df = reduce(df, migrated_cgIDs)

        return df

    def chunkIt(seq, num):
        """
        Split the feature sets into chunks
        """
        avg = len(seq) / float(num)
        out = []
        last = 0.0

        while last < len(seq):
            out.append(seq[int(last) : int(last + avg)])
            last += avg

        return out

    def CT(self, file, molecule):
        if molecule == "RNA":
            r = re.compile("TCGA-[a-zA-Z\d]{2,6}")
            m = r.search(file)
            if m:
                cancer_type = str(m.group(0))
            else:
                r = re.compile("TARGET-[a-zA-Z\d]{2,6}")
                m = r.search(file)
                cancer_type = str(m.group(0))

        if molecule == "DNA":
            """
            Str extract the Cancer type for the
            """
            r = re.compile("(?<=level_)[a-zA-Z\d]{3,8}")
            m = r.search(file)
            if m:
                cancer_type = str(m.group(0))

        if molecule == "methylation":
            r = re.compile("TCGA-[a-zA-Z\d]{2,4}")
            m = r.search(file)

            if m:
                cancer_type = str(m.group(0))
        return cancer_type

    def labeler(self, df, cancer_type, molecule):
        if molecule == "DNA":
            labels = []
            index = df.index
            sample_list = list(index)
            for sample in sample_list:
                # extract the last 3 chars in the sample ID column
                r = re.compile("(?<=-)\d\d")
                m = r.search(str(sample))
                if m:
                    if "01" in str(m.group(0)) or "11" in str(m.group(0)):
                        labels.append("Tumor")
                    else:
                        labels.append("Normal")
                else:
                    labels.append("remove")
            # Now the labels column to df as the labels column
            df["cancerType"] = cancer_type
            # add the cancertype column
            df["labels"] = labels

        if molecule == "RNA":
            labels = []
            index = df.index
            sample_list = list(index)
            for sample in sample_list:
                # extract the last 3 chars in the sample ID column
                r = re.compile("(?<=-)\d\d")
                m = r.search(str(sample))
                if m:
                    if "01" in str(m.group(0)) or "11" in str(m.group(0)):
                        labels.append("Tumor")
                    else:
                        labels.append("Normal")
                else:
                    labels.append("remove")
            # Now the labels column to df as the labels column
            df["cancerType"] = cancer_type
            # add the cancertype column
            df["labels"] = labels

        if molecule == "methylation":
            labels = []
            index = df.index
            sample_list = list(index)
            for sample in sample_list:
                # extract the last 3 chars in the sample ID column
                r = re.compile("(?<=-)\d\dA")
                m = r.search(str(sample))
                if m:
                    if "01" in str(m.group(0)) or "11" in str(m.group(0)):
                        labels.append("Tumor")
                    else:
                        labels.append("Normal")
                else:
                    labels.append("remove")
            # Now the labels column to df as the labels column
            df["cancerType"] = cancer_type
            # add the cancertype column
            df["labels"] = labels
        return df

    def build_input(self, path, start, end, molecule):
        block = []
        for file in glob.glob(path):
            dat = test.read_cut(file, start, end, molecule)
            cancer_type = test.CT(file, molecule)
            df = test.labeler(dat, cancer_type, molecule)
            block.append(df)
        df_block = pd.concat(block, ignore_index=True)
        df_block = df_block[df_block.labels != "remove"]

        return df_block

    def feature_selection(self, X, molecule):
        pass

    def grade_all(self, path, molecule):
        pass

    def join(self, path, molecule):
        pass

    def subset(self, path, features, molecule):
        pass

    def encode(self, y):
        pass

    def data_preprocessing(self, X, molecule):
        pass

    def balance(self, X, y, molecule):
        pass

    def evaluate(self, X, y, model, molecule):
        pass


class Classifier:
    def build_calssifier(self):
        pass

    def split(self):
        pass

    def fit_clf(self):
        pass

    def visualize(self):
        pass


# In[104]:


# https://scikit-learn.org/stable/modules/ensemble.html


# In[105]:


test = molecule_preprocessing(path="/home/jovyan/CSBL_shared/RNASeq/TCGA/counts/*.csv")


# 
# The random forest regression seems not to be the optimal option
# the two data type feature selection. We will need to regroup here for better perfromance.
# Further the sequential NN may not be best suited for the applicaitons here.
# If we intend to cap the DNNs with a voting classifier is may be pertanent to use the 
# functional NN that is more fluid with TF
# 
# 
# 
# 
# First three approaches:
# 1. Filter based: 
#     We specify some metric and based on that filter features. 
#     An example of such a metric could be correlation/chi-square.
# 2. Wrapper-based: 
#     Wrapper methods consider the selection of a set of 
#     features as a search problem. Example: Recursive 
#     Feature Elimination
# 3. Embedded: Embedded methods use algorithms that have 
#     built-in feature selection methods. For instance, 
#     Lasso and RF have their own feature selection methods.
# 

# In[113]:


# Correlation
class feature_selection:
    def __init__(self, path, X):
        self.path = path
        self.X = X

        pass

    def encode_y(y):
        y = np.array(y)
        y = y.reshape(-1, 1)
        encoder = OneHotEncoder(sparse=False)
        encoder.fit(y)
        y_encoded = encoder.transform(y)

        return y_encoded
    
    def dummies(X):
        if "Composite Element REF" in X.columns:
            X = X.drop("Composite Element REF", axis=1)
        if "samples" in X.columns:
            X = X.drop("samples", axis=1)
        X_processed = pd.get_dummies(X, columns=["labels"])

        return X_processed

    def cor_selector(X, num_feats):
        X = feature_selection.dummies(X)
        y = X["cancerType"]
        X = X.drop("cancerType", axis=1)
        cor_list = []
        feature_name = X.columns.tolist()
        # calculate the correlation with y for each feature
        for i in X.columns.tolist():
            cor = np.corrcoef(X[i], y)[0, 1]
            cor_list.append(cor)
        # replace NaN with 0
        cor_list = [0 if np.isnan(i) else i for i in cor_list]
        # feature name
        cor_feature = X.iloc[
            :, np.argsort(np.abs(cor_list))[-num_feats:]
        ].columns.tolist()
        # feature selection? 0 for not select, 1 for select
        cor_support = [True if i in cor_feature else False for i in feature_name]
        return cor_support, feature_name

    # cor_support, feature_name = cor_selector(X, num_feats=5000)

    def chi_selector(X, num_feats):
        X = feature_selection.dummies(X)
        y = X["cancerType"]
        X = X.drop("cancerType", axis=1)
        X_norm = MinMaxScaler().fit_transform(X)
        chi_selector = SelectKBest(chi2, k=num_feats)
        chi_selector.fit(X_norm, y)
        chi_support = chi_selector.get_support()
        chi_feature = X.loc[:, chi_support].columns.tolist()

        return chi_support
    
    def rfR(X, num_feats):
        X = feature_selection.dummies(X)
        y = X["cancerType"]
        X = X.drop("cancerType", axis=1)
        y_encoded = feature_selection.encode_y(y)
        X_norm = MinMaxScaler().fit_transform(X)
        X_train,y_train,X_test,y_test = train_test_split(X_norm,test_size=0.3)
        sel = SelectFromModel(estimator =RandomForestRegressor(n_estimators = 100,n_jobs=20))
        sel.fit(X_train, y_train)
        rfR_support = sel.get_support()
        rfe_feature = X_train.columns[(sel.get_support())]
        
        return rfR_support

    # chi_support = chi_selector(X, num_feats=5000)

    def logR(X, num_features):
        X = feature_selection.dummies(X)
        y = X["cancerType"]
        X = X.drop("cancerType", axis=1)
        y_encoded = feature_selection.encode_y(y)
        X_norm = MinMaxScaler().fit_transform(X)
        rfe_selector = RFE(
            estimator=LogisticRegression(),
            n_features_to_select=num_feats,
            step=10,
            verbose=5,
        )
        rfe_selector.fit(X_norm, y)
        rfe_support = rfe_selector.get_support()
        rfe_feature = X.loc[:, rfe_support].columns.tolist()

        return rfe_support

    # rfe_support = logR(X, num_feats=5000)

    def lassoR(X, num_feats):
        X = feature_selection.dummies(X)
        y = X["cancerType"]
        X = X.drop("cancerType", axis=1)
        y_encoded = feature_selection.encode_y(y)
        X_norm = MinMaxScaler().fit_transform(X)
        embeded_lr_selector = SelectFromModel(
            LogisticRegression(penalty="l1"), max_features=num_feats
        )
        embeded_lr_selector.fit(X_norm, y)
        embeded_lr_support = embeded_lr_selector.get_support()
        embeded_lr_feature = X.loc[:, embeded_lr_support].columns.tolist()
        return embeded_lr_support, embeded_lr_feature

    # embeded_lr_support = lassoR(X, num_feats=5000)

    def rfC(X, num_feats):
        X = feature_selection.dummies(X)
        y = X["cancerType"]
        X = X.drop("cancerType", axis=1)
        y_encoded = feature_selection.encode_y(y)
        X_norm = MinMaxScaler().fit_transform(X)
        embeded_rf_selector = SelectFromModel(
            RandomForestClassifier(n_estimators=100), max_features=num_feats
        )
        embeded_rf_selector.fit(X, y)
        embeded_rf_support = embeded_rf_selector.get_support()
        embeded_rf_feature = X.loc[:, embeded_rf_support].columns.tolist()

        return embeded_rf_support

    # embeded_rf_support = rfC(X, num_feats=5000)

    def LGBMC_selector(X, num_feats):
        X = feature_selection.dummies(X)
        y = X["cancerType"]
        X = X.drop("cancerType", axis=1)
        y_encoded = feature_selection.encode_y(y)
        X_norm = MinMaxScaler().fit_transform(X)
        lgbc = LGBMClassifier(
            n_estimators=500,
            learning_rate=0.05,
            num_leaves=32,
            colsample_bytree=0.2,
            reg_alpha=3,
            reg_lambda=1,
            min_split_gain=0.01,
            min_child_weight=40,
        )

        embeded_lgb_selector = SelectFromModel(lgbc, max_features=num_feats)
        embeded_lgb_selector.fit(X_norm, y)

        embeded_lgb_support = embeded_lgb_selector.get_support()
        embeded_lgb_feature = X.loc[:, embeded_lgb_support].columns.tolist()

        return embeded_lgb_support

    # embeded_lgb_support = rfC(X, num_feats=5000)

    def cross_validate_feature_selection(
        feature_name,
        cor_support,
        chi_support,
        rfe_support,
        embedded_lr_support,
        embedded_rf_support,
        embedded_lgb_support,
    ):
        df = pd.DataFrame(
            {
                "Feature": feature_name,
                "Pearson": cor_support,
                "Chi-2": chi_support,
                "RFE": rfe_support,
                "Logistics": embeded_lr_support,
                "RandomForestClassifier": embeded_rf_support,
                "RandomForstRegression "
                "LightGBM": embeded_lgb_support,
            }
        )
        # count the selected times for each feature
        df["Total"] = np.sum(df, axis=1)
        df = df.sort_values(["Total", "Feature"], ascending=False)
        df.index = range(1, len(df) + 1)

        return df
    
    def filter_features(X, feature_selection_df):
        # select 5000 features from the X data frames
        # return these for modeling the data with the DNNs
        pass


# In[107]:


# feature_selection_df = cross_val_score(
#    feature_name,
#    cor_support,
#    chi_support,
#    rfe_support,
#    embedded_lr_support,
#    embedded_rf_support,
#    embedded_lgb_support,
# )

# feature_selection_df.to_csv("window_" + str(feature_name[0]) + "_" + str(feature_name[-1]) + ".csv")


# In[108]:


print(dir(molecule_preprocessing))


# In[109]:


df = test.read_cut(
    file="/home/jovyan/CSBL_shared/RNASeq/TCGA/counts/TCGA-SARC.counts.csv",
    start=0,
    end=100,
    molecule="RNA",
)
# cancer_type = test.CT(file = "/home/jovyan/CSBL_shared/RNASeq/TCGA/counts/TCGA-SARC.counts.csv",molecule = "RNA")
# df = test.labeler(df, cancer_type, molecule ="RNA")
df = test.build_input(
    path="/home/jovyan/CSBL_shared/RNASeq/TCGA/counts/*.csv",
    start=0,
    end=604,
    molecule="RNA",
)


# In[114]:


df


# In[ ]:


# get feature importances and take top 50

