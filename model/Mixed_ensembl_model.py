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

    def build_input(self, path, start, end):
        block = []
        for file in glob.glob(path):
            dat = read_cut(file, start, end)
            cancer_type = CT(file)
            df = labeler(dat, cancer_type)
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
