# import the necessary libraries
import pandas as pd
import numpy as np
import re
import os

# read in file as a pandas data frame
directory = "/Volumes/easy store/Cornell/Methylation_project/Expanded/"
os.chdir(directory)
for filename in os.listdir(directory):
    if(filename.endswith(".methylation450.tsv")):
        print("Reading:" + filename)

        # Read in the file as a pandas data frame
        dat = pd.read_table(filename)
        print(dat.shape)

        # transpose the data frame:
        dat_t = dat.T
        print(dat_t.shape)

        # extract the label value we want to rep in the column named labels
        txt = filename
        x = re.search("\w{4}-\w{2,5}", txt)
        print(x.group())

        # make an array of the length of the rows, fill it with one value
        print(dat_t[0].count())
        val = np.array([x.group()])
        col = np.repeat(val, dat_t[0].count(), axis=0)

        # add the array for as labels column of  the transposed data frame:
        dat_t['labels'] = col
        print("labels attached")
        # write the file
        print("Writing labeled and transposed: " + x.group())
        dat_t.to_csv(str("/Volumes/easy store/Cornell/Methylation_project/Expanded/" + x.group() + ".csv"), index=False)

'''
# one at a time for testing steps
dat = pd.read_table("/Volumes/easy store/Cornell/Methylation_project/Expanded/TCGA-ACC.methylation450.tsv")
print(dat.shape)

# transpose the data frame:

dat_t = dat.T
print(dat_t.shape)

# extract the label value we want to rep in the column named labels
txt = "/Volumes/easy store/Cornell/Methylation_project/Expanded/TCGA-ACC.methylation450.tsv"
x = re.search("\w{4}-\w{2,5}", txt)
print(x.group())

# make an array of the length of the rows, fill it with one value
print(dat_t[0].count())
val = np.array([x.group()])
col = np.repeat(val, dat_t[0].count(), axis=0)

# add the array for as labels column of  the transposed data frame:
dat_t['labels'] = col

# write the file
#compression_opts = dict(method='zip', archive_name=str("/Volumes/easy store/Cornell/Methylation_project/Expanded/" +
#                                                       x.group() + ".csv"))"

dat_t.to_csv(str("/Volumes/easy store/Cornell/Methylation_project/Expanded/" + x.group() + ".csv"), index=False)

# now put it all together
'''





