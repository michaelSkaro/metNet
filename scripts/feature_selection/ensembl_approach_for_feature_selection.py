# add an ensemble method for the RNA and deal with the DNA
"""
The random forest regression seems not to be the optimal option
the two data type feature selection. We will need to regroup here for better perfromance.
Further the sequential NN may not be best suited for the applicaitons here.
If we intend to cap the DNNs with a voting classifier is may be pertanent to use the 
functional NN that is more fluid with TF
"""


"""
First three approaches:
1. Filter based: 
    We specify some metric and based on that filter features. 
    An example of such a metric could be correlation/chi-square.
2. Wrapper-based: 
    Wrapper methods consider the selection of a set of 
    features as a search problem. Example: Recursive 
    Feature Elimination
3. Embedded: Embedded methods use algorithms that have 
    built-in feature selection methods. For instance, 
    Lasso and RF have their own feature selection methods.
"""


from lightgbm import LGBMClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import SelectFromModel, SelectKBest, chi2
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import MinMaxScaler


# Correlation
def cor_selector(X, y, num_feats):
    y = y = X["cancerType"]
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
    cor_feature = X.iloc[:, np.argsort(np.abs(cor_list))[-num_feats:]].columns.tolist()
    # feature selection? 0 for not select, 1 for select
    cor_support = [True if i in cor_feature else False for i in feature_name]
    return cor_support, cor_feature, feature_name


# cor_support, cor_feature, feature_name = cor_selector(X, num_feats=5000)


def chi_selector(X, num_feats):
    y = y = X["cancerType"]
    X = X.drop("cancerType", axis=1)
    X_norm = MinMaxScaler().fit_transform(X)
    X_filtered = SelectKBest(chi2, k=num_feats).fit_transform(X_norm, y)
    return X_filtered


# X_filtered = chi_selector(X, num_feats =1000)


def logR(X, num_features):
    y = y = X["cancerType"]
    X = X.drop("cancerType", axis=1)
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

    return rfe_support, rfe_feature


# rfe_support, rfe_feature = logR(X, num_feats=5000)


def lassoR(X, num_feats):
    y = y = X["cancerType"]
    X = X.drop("cancerType", axis=1)
    X_norm = MinMaxScaler().fit_transform(X)
    embeded_lr_selector = SelectFromModel(
        LogisticRegression(penalty="l1"), max_features=num_feats
    )
    embeded_lr_selector.fit(X_norm, y)
    embeded_lr_support = embeded_lr_selector.get_support()
    embeded_lr_feature = X.loc[:, embeded_lr_support].columns.tolist()
    return embeded_lr_support, embeded_lr_feature


# embeded_lr_support, embeded_lr_feature = lassoR(X, num_feats=5000)


def rfC(X, num_feats):
    y = y = X["cancerType"]
    X = X.drop("cancerType", axis=1)
    X_norm = MinMaxScaler().fit_transform(X)
    embeded_rf_selector = SelectFromModel(
        RandomForestClassifier(n_estimators=100), max_features=num_feats
    )
    embeded_rf_selector.fit(X, y)
    embeded_rf_support = embeded_rf_selector.get_support()
    embeded_rf_feature = X.loc[:, embeded_rf_support].columns.tolist()

    return embeded_rf_support, embeded_rf_feature


# embeded_rf_support, embeded_rf_feature = rfC(X, num_feats=5000)


def LGBMC_selector(X, num):
    y = y = X["cancerType"]
    X = X.drop("cancerType", axis=1)
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

    return embeded_lgb_support, embeded_lgb_feature


# embeded_lgb_support, embeded_lgb_feature = rfC(X, num_feats=5000)


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
            "Random Forest": embeded_rf_support,
            "LightGBM": embeded_lgb_support,
        }
    )
    # count the selected times for each feature
    df["Total"] = np.sum(df, axis=1)
    df = df.sort_values(["Total", "Feature"], ascending=False)
    df.index = range(1, len(df) + 1)

    return df


# feature_selection_df = cross_val_score(feature_namer, cor_support, chi_support, rfe_support, embedded_lr_support, embedded_rf_support, embedded_lgb_support)
