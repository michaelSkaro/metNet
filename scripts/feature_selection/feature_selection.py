# Correlation
class feature_selection:
    def __init__(self, path, X):
        self.path = path
        self.X = X

        pass

    def cor_selector(X, num_feats):
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
        cor_feature = X.iloc[
            :, np.argsort(np.abs(cor_list))[-num_feats:]
        ].columns.tolist()
        # feature selection? 0 for not select, 1 for select
        cor_support = [True if i in cor_feature else False for i in feature_name]
        return cor_support, feature_name

    cor_support, feature_name = cor_selector(X, num_feats=5000)

    def chi_selector(X, num_feats):
        y = y = X["cancerType"]
        X = X.drop("cancerType", axis=1)
        X_norm = MinMaxScaler().fit_transform(X)
        chi_selector = SelectKBest(chi2, k=num_feats)
        chi_selector.fit(X_norm, y)
        chi_support = chi_selector.get_support()
        chi_feature = X.loc[:, chi_support].columns.tolist()

        return chi_support

    chi_support = chi_selector(X, num_feats=1000)

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

        return rfe_support

    rfe_support = logR(X, num_feats=5000)

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

    embeded_lr_support = lassoR(X, num_feats=5000)

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

        return embeded_rf_support

    embeded_rf_support = rfC(X, num_feats=5000)

    def LGBMC_selector(X, num_feats):
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

        return embeded_lgb_support

    embeded_lgb_support = rfC(X, num_feats=5000)

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

    feature_selection_df = cross_val_score(
        feature_name,
        cor_support,
        chi_support,
        rfe_support,
        embedded_lr_support,
        embedded_rf_support,
        embedded_lgb_support,
    )
