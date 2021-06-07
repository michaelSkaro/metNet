#!/usr/bin/env python
# coding: utf-8

# In[469]:


# New approach lets attack building the shell.
import glob
import io
import math
import os
import re
import time
import warnings
from collections import Counter

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tensorflow as tf
from imblearn.over_sampling import SMOTE
from lightgbm import LGBMClassifier
from matplotlib import pyplot
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.feature_selection import RFE, SelectFromModel, SelectKBest, chi2
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import explained_variance_score
from sklearn.model_selection import cross_val_score, train_test_split
from sklearn.preprocessing import (
    LabelEncoder,
    MinMaxScaler,
    OneHotEncoder,
    OrdinalEncoder,
)
from tensorflow import keras
from tensorflow.keras import Sequential, backend, layers
from tensorflow.keras.layers import Dense, Dropout


# In[605]:


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
            df = df.drop("sample", axis=1)

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
            dat = mp.read_cut(file, start, end, molecule)
            cancer_type = mp.CT(file, molecule)
            df = mp.labeler(dat, cancer_type, molecule)
            block.append(df)
        df_block = pd.concat(block, ignore_index=True)
        df_block = df_block[df_block.labels != "remove"]

        return df_block

    def join(self, path, molecule):

        if molecule == "RNA" or "DNA":
            features = []
            for file in glob.glob(path):
                dat = pd.read_csv(file)
                features.append(dat)
                # concatnenate them and sort the feature ranks
            df = pd.concat(features, ignore_index=True)
            df.drop_duplicates(subset=["Feature"])
            df.sort_values(by="Total", ascending=False)

            features = df.Feature
            features = features.drop_duplicates()

        if molecule == "methylation":
            features = []
            for file in glob.glob(path):
                dat = pd.read_csv(file)
                features.append(dat)
                # concatnenate them and sort the feature ranks
            df = pd.concat(features, ignore_index=True)
            df.drop_duplicates(subset=["cgIDs"])
            df.sort_values(by="importance_scores", ascending=False)
            # get cgIDs
            features = df.cgIDs
            features = features.drop_duplicates()

        return features

    def subset(self, path, features, molecule):
        selected = []
        tic_all = time.perf_counter()
        for file in glob.glob(path):
            print(file)
            df = pd.read_table(file, delimiter=",")
            df = df[df["Ensembl"].isin(features)]
            df = df.T
            names = df.iloc[0]
            cancer_type = molecule_preprocessing.CT(self, file, molecule)
            df = molecule_preprocessing.labeler(self, df, cancer_type, molecule)
            print(df.shape)
            selected.append(df)
        df_block = pd.concat(selected, ignore_index=True)
        df_block = df_block.fillna(0)
        df_block = df_block[df_block.labels != "remove"]
        print(df_block.shape)
        column_indices = range(4707)
        new_names = names
        old_names = df_block.columns[column_indices]
        df_block.rename(columns=dict(zip(old_names, new_names)), inplace=True)
        toc_all = time.perf_counter()
        print(f"Read in features in {toc_all - tic_all:0.4f} seconds")
        print(features.shape)
        return df_block

    def DNN_preprocessing(self, X):

        # if the file has rowIDs, we can drop these
        if "Composite Element REF" in X.columns:
            X = X.drop("Composite Element REF", axis=1)
        if "Unnamed: 0" in X.columns:
            X = X.drop("Unnamed: 0", axis=1)
        if "sample" in X.columns:
            X = X.drop("sample", axis=1)
        # subset out the rows that have Target
        X = X[X.cancerType != "TARGET-WT"]
        X = X[X.cancerType != "TARGET-RT"]
        X = X[X.cancerType != "TARGET-NBL"]
        X = X[X.cancerType != "TARGET-AML"]
        # combined READ and COAD
        X = X.replace({"cancerType": "READ"}, "COADREAD")
        X = X.replace({"cancerType": "COAD"}, "COADREAD")
        # substring the TCGA- out of the column
        # get rid of TCGA- in cancer type
        X["cancerType"] = X["cancerType"].str[5:]
        y = X["cancerType"]
        X = X.drop("cancerType", axis=1)
        X = pd.get_dummies(X, columns=["labels"])
        return X, y

    def balance(self, X, y):
        y_encoded = LabelEncoder().fit_transform(y)
        counter = Counter(y_encoded)
        for k, v in counter.items():
            per = v / len(y) * 100
            print("Class=%d, n=%d (%.3f%%)" % (k, v, per))
        # plot the distribution
        pyplot.bar(counter.keys(), counter.values())
        pyplot.show()

        oversample = SMOTE()
        X, y = oversample.fit_resample(X, y_encoded)
        # summarize distribution
        counter = Counter(y)
        for k, v in counter.items():
            per = v / len(y) * 100
            print("Class=%d, n=%d (%.3f%%)" % (k, v, per))
        # plot the distribution
        pyplot.bar(counter.keys(), counter.values())
        pyplot.show()
        for column in X:
            X[column] = X[column].apply(np.ceil)
        X = X.astype(int)

        return X, y


# In[606]:


# Correlation
class feature_selection:
    def __init__(self, X):
        self.X = X

        pass

    def encode_y(y):
        y = np.array(y)
        y = y.reshape(-1, 1)
        encoder = OrdinalEncoder()
        # encoder = OneHotEncoder(sparse=False)
        encoder.fit(y)
        y_encoded = encoder.transform(y)

        return y_encoded

    def prep(X):
        y = X["cancerType"]
        X = X.drop("cancerType", axis=1)
        y = feature_selection.encode_y(y)
        X = feature_selection.dummies(X)
        feature_list = list(X.columns)
        # X = MinMaxScaler().fit_transform(X)

        return X, y, feature_list

    def dummies(X):
        if "Composite Element REF" in X.columns:
            X = X.drop("Composite Element REF", axis=1)
        if "samples" in X.columns:
            X = X.drop("samples", axis=1)
        X = pd.get_dummies(X, columns=["labels"])

        return X

    def chi_selector(X, num_feats):
        X, y, feature_list = feature_selection.prep(X)
        y = feature_selection.encode_y(y)
        chi_selector = SelectKBest(chi2, k=num_feats)
        chi_selector.fit(X, y.ravel())
        chi_support = chi_selector.get_support()
        # chi_feature = X.loc[:, chi_support].columns.tolist()

        return chi_support, feature_list

    def rfR(X, num_feats):
        X, y, feature_list = feature_selection.prep(X)
        y = feature_selection.encode_y(y)
        embeded_rf_selector = SelectFromModel(
            RandomForestRegressor(n_estimators=100, n_jobs=20), max_features=num_feats
        )
        embeded_rf_selector.fit(X, y.ravel())
        embeded_rf_support = embeded_rf_selector.get_support()
        # embeded_rf_feature = X.loc[:, embeded_rf_support].columns.tolist()

        return embeded_rf_support

    # chi_support = chi_selector(X, num_feats=5000)

    def logR(X, num_feats):
        X, y, feature_list = feature_selection.prep(X)
        y = feature_selection.encode_y(y)
        rfe_selector = RFE(
            estimator=LogisticRegression(n_jobs=20),
            n_features_to_select=num_feats,
            step=10,
            verbose=5,
        )
        rfe_selector.fit(X, y.ravel())
        rfe_support = rfe_selector.get_support()
        # rfe_feature = X.loc[:, rfe_support].columns.tolist()

        return rfe_support

    # rfe_support = logR(X, num_feats=5000)

    def lassoR(X, num_feats):
        X, y, feature_list = feature_selection.prep(X)
        embeded_lr_selector = SelectFromModel(
            LogisticRegression(penalty="l2", n_jobs=20), max_features=num_feats
        )
        embeded_lr_selector.fit(X, y.ravel())
        embeded_lr_support = embeded_lr_selector.get_support()
        # embeded_lr_feature = X.loc[:, embeded_lr_support].columns.tolist()
        return embeded_lr_support

    # embeded_lr_support = lassoR(X, num_feats=5000)

    def rfC(X, num_feats):
        X, y, feature_list = feature_selection.prep(X)
        embeded_rf_selector = SelectFromModel(
            RandomForestClassifier(n_estimators=100, n_jobs=20), max_features=num_feats
        )
        embeded_rf_selector.fit(X, y.ravel())
        embeded_rf_support = embeded_rf_selector.get_support()
        # embeded_rf_feature = X.loc[:, embeded_rf_support].columns.tolist()

        return embeded_rf_support

    # embeded_rf_support = rfC(X, num_feats=5000)

    def LGBMC_selector(X, num_feats):
        X, y = feature_selection.prep(X)
        lgbc = LGBMClassifier(
            n_estimators=500,
            learning_rate=0.05,
            num_leaves=32,
            colsample_bytree=0.2,
            reg_alpha=3,
            reg_lambda=1,
            min_split_gain=0.01,
            min_child_weight=40,
            n_jobs=20,
        )

        embeded_lgb_selector = SelectFromModel(lgbc, max_features=num_feats)
        embeded_lgb_selector.fit(X, y.ravel())

        embeded_lgb_support = embeded_lgb_selector.get_support()
        # embeded_lgb_feature = X.loc[:, embeded_lgb_support].columns.tolist()

        return embeded_lgb_support

    # embeded_lgb_support = rfC(X, num_feats=5000)

    def cross_validate_feature_selection(
        feature_list,
        chi_support,
        rfe_support,
        embeded_lr_support,
        embeded_rfC_support,
        embeded_rfR_support,
        # embeded_lgb_support,
    ):
        df = pd.DataFrame(
            {
                "Feature": feature_list,
                "Chi-2": chi_support,
                "RFE": rfe_support,
                "Logistics": embeded_lr_support,
                "RandomForestClassifier": embeded_rfC_support,
                "RandomForstRegression": embeded_rfR_support,
                # "LightGBM": embeded_lgb_support, # try to fix this later?
            }
        )
        # count the selected times for each feature
        df["Total"] = np.sum(df, axis=1)
        df = df.sort_values(["Total", "Feature"], ascending=False)
        df.index = range(1, len(df) + 1)

        return df

    def grade_features(X):
        # cor_support, feature_list = feature_selection.cor_selector(
        #    X, num_feats=50)
        chi_support, feature_list = feature_selection.chi_selector(X, num_feats=50)
        rfe_support = feature_selection.lassoR(X, num_feats=50)
        embeded_lr_support = feature_selection.logR(X, num_feats=50)
        embeded_rfC_support = feature_selection.rfC(X, num_feats=50)
        embeded_rfR_support = feature_selection.rfR(X, num_feats=50)
        # embeded_lgb_support = feature_selection.LGBMC_selector(X, num_feats=50)

        CV = feature_selection.cross_validate_feature_selection(
            feature_list,
            # cor_support,
            chi_support,
            rfe_support,
            embeded_lr_support,
            embeded_rfC_support,
            embeded_rfR_support,
            # embeded_lgb_support,
        )
        CV = CV[1:50]

        return CV


# In[ ]:


def grade_all(path, molecule):
    if molecule == "RNA":
        windows = molecule_preprocessing.chunkIt(seq=range(60486), num=100)

        for i in windows:
            start = list(i)[0]
            end = list(i)[-1]
            df_block = molecule_preprocessing.build_input(
                mp, path=path, start=start, end=end, molecule="RNA"
            )
            CV = feature_selection.grade_features(df_block)
            CV.to_csv("Ensembl_ID_window_" + str(start) + "_" + str(end) + ".csv")

    if molecule == "DNA":
        windows = molecule_preprocessing.chunkIt(seq=range(40545), num=100)

        for i in windows:
            start = list(i)[0]
            end = list(i)[-1]
            df_block = molecule_preprocessing.build_input(
                mp, path=path, start=start, end=end, molecule="DNA"
            )
            CV = feature_selection.grade_features(df_block)
            CV.to_csv("Gene_ID_window_" + str(start) + "_" + str(end) + ".csv")

    if molecule == "methylation":
        windows = molecule_preprocessing.chunkIt(seq=range(485577), num=100)

        for i in windows:
            start = list(i)[0]
            end = list(i)[-1]
            df_block = molecule_preprocessing.build_input(
                mp, path=path, start=start, end=end, molecule="methylation"
            )
            CV = feature_selection.grade_features(df_block)
            CV.to_csv("Methylation_cgID_window_" + str(start) + "_" + str(end) + ".csv")

    pass


df = test.build_input(
    path="/home/jovyan/CSBL_shared/RNASeq/TCGA/counts/*.csv",
    start=0,
    end=604,
    molecule="RNA",
)
CV = feature_selection.grade_features(df)
mp = molecule_preprocessing(path="/home/jovyan/CSBL_shared/RNASeq/TCGA/counts/*.csv")
fs = feature_selection(mp)
df = grade_all(path="/home/jovyan/CSBL_shared/RNASeq/TCGA/counts/*.csv", molecule="RNA")
mp = molecule_preprocessing(
    path="/home/jovyan/storage/Machine_Learning/RNA_Selected_IDs_ensemble/"
)
features = molecule_preprocessing.join(
    self=mp,
    path="/home/jovyan/storage/Machine_Learning/RNA_Selected_IDs_ensemble/*.csv",
    molecule="RNA",
)
df = molecule_preprocessing.subset(
    self=mp,
    path="/home/jovyan/CSBL_shared/RNASeq/TCGA/counts/*.csv",
    features=features,
    molecule="RNA",
)
df.to_csv("Selected_features_all_RNA_data_TCGA_6_4_21.csv")
X, y = molecule_preprocessing.DNN_preprocessing(mp, df)
X, y = molecule_preprocessing.balance(mp, X, y)


# In[478]:


class Classifier:
    def build_calssifier(self):
        pass

    def split(self):
        pass

    def fit_clf(self):
        pass

    def visualize(self):
        pass


# In[537]:


df = molecule_preprocessing.subset(
    self=mp,
    path="/home/jovyan/CSBL_shared/RNASeq/TCGA/counts/*.csv",
    features=features,
    molecule="RNA",
)
df.to_csv("Selected_features_all_RNA_data_TCGA_6_4_21.csv")
X, y = molecule_preprocessing.DNN_preprocessing(mp, df)
X, y = molecule_preprocessing.balance(mp, X, y)


# In[541]:


from numpy import loadtxt
from sklearn.metrics import accuracy_score
from sklearn.model_selection import train_test_split
from xgboost import XGBClassifier

X.astype("float16").dtypes
seed = 19  # love you bubba
test_size = 0.33
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=test_size, random_state=seed
)


# In[542]:


model = XGBClassifier()
model.fit(X_train, y_train.ravel())


# In[545]:


X.astype("float16").dtypes


# In[549]:


X = X.astype(int)
X.dtypes


# In[551]:


seed = 19  # love you bubba
test_size = 0.33
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=test_size, random_state=seed
)
model = XGBClassifier(n_jobs=30)
model.fit(X_train, y_train.ravel())


# In[552]:


# make predictions for test data
y_pred = model.predict(X_test)
predictions = [round(value) for value in y_pred]
# evaluate predictions
accuracy = accuracy_score(y_test, predictions)
print("Accuracy: %.2f%%" % (accuracy * 100.0))


# In[ ]:


"the default evaluation metric used with the objective 'multi:softprob' was changed from 'merror' to 'mlogloss'""


# In[564]:


# Plot non-normalized confusion matrix
from sklearn.metrics import (
    classification_report,
    plot_confusion_matrix,
    plot_det_curve,
    plot_roc_curve,
)

titles_options = [
    ("Confusion matrix, without normalization", None),
    ("Normalized confusion matrix", "true"),
]
for title, normalize in titles_options:
    disp = plot_confusion_matrix(
        model, X_test, y_test, cmap=plt.cm.Blues, normalize=normalize
    )
    disp.ax_.set_title(title)

    print(title)
    print(disp.confusion_matrix)

plt.show()


# In[566]:


# Do DNA


# In[ ]:


mp = molecule_preprocessing(
    path="/home/jovyan/storage/UCSC_Xena/somatic_variation/UCSC_TCGA_mutations/*.gz"
)
#df = molecule_preprocessing.build_input(self =mp, path="/home/jovyan/storage/UCSC_Xena/somatic_variation/UCSC_TCGA_mutations/*.gz",start=0,end=405,molecule="DNA")
#CV = feature_selection.grade_features(df)
fs = feature_selection(mp)
df = grade_all(path="/home/jovyan/storage/UCSC_Xena/somatic_variation/UCSC_TCGA_mutations/*.gz", molecule="DNA")


# In[ ]:





# In[ ]:


mp = molecule_preprocessing(
    path="/home/jovyan/storage/Machine_Learning/RNA_Selected_IDs_ensemble/"
)
features = molecule_preprocessing.join(
    self=mp,
    path="/home/jovyan/storage/Machine_Learning/RNA_Selected_IDs_ensemble/*.csv",
    molecule="RNA",
)
df = molecule_preprocessing.subset(
    self=mp,
    path="/home/jovyan/CSBL_shared/RNASeq/TCGA/counts/*.csv",
    features=features,
    molecule="RNA",
)
df.to_csv("Selected_features_all_RNA_data_TCGA_6_4_21.csv")
X, y = molecule_preprocessing.DNN_preprocessing(mp, df)
X, y = molecule_preprocessing.balance(mp, X, y)

