import numpy as np
#import statsmodels.api as sm

import pandas as pd
 
train_df = pd.read_csv('energies.csv')

example_series =  pd.Series([1,5,10,30,50,30,15,40,45])
print(train_df['NB'].median())
print(train_df.cov())

print(train_df.corr())



print(train_df.columns)
#train_df['RMSD'].describe()
#dfols=train_df[train_df['RMSD']]






