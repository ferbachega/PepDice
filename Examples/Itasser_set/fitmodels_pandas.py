import numpy as np
#import statsmodels.api as sm

import pandas as pd
 
train_df = pd.read_csv('energies_gyration_AB_energy_new.csv')

example_series =  pd.Series([1,5,10,30,50,30,15,40,45])
print(train_df['NB'].median())
print(train_df.cov())

print(train_df.corr())

allkeys = []
keys    = [] 
wkeys   = []

for key in train_df:
    print key 
    if key == 'PDB' or key == 'decoy' or key =='RMSD' or key == 'CONTACT' or key =='R_GYRATION' or key == 'AB_ENERGY' :
        pass
    else:
        train_df[key+'/size'] = train_df[key]/train_df['SIZE']
        allkeys.append(key+'/size')
        allkeys.append(key)
        
        keys.append(key)
        wkeys.append(key+'/size')


print(train_df.columns)
#train_df['RMSD'].describe()
#dfols=train_df[train_df['RMSD']]



pdbs = [
        'IT1af7__' ,
        'IT1ah9_'  ,
        'IT1aoy_'  ,
        'IT1b4bA'  ,
        'IT1b72A'  ,
        'IT1bm8_'  ,
        'IT1dcjA_' ,
        'IT1dtjA_' ,
        'IT1egxA'  ,
        'IT1fo5A'  ,
        'IT1g1cA'  ,
        'IT1gjxA'  ,
        'IT1gpt_'  ,
        'IT1gyvA'  ,
        'IT1itpA'  ,
        'IT1kjs_'  ,
        'IT1kviA'  ,
        'IT1mkyA3' ,
        'IT1mla_2' ,
        'IT1n0uA4' ,
        'IT1ne3A'  ,
        'IT1npsA'  ,
        'IT1o2fB_' ,
        'IT1of9A'  ,
        'IT1r69_'  ,
        'IT1shfA'  ,
        'IT1sro_'  ,
        'IT1tfi_'  ,
        'IT1tif_'  ,
        'IT1tig_'  ,
        'IT1vcc_'  ,
        'IT2cr7A'  ,
        'IT2f3nA'  ,
        'IT2pcy_'  ,
        'IT2reb_2' ,
        'IT256bA'  ,
        ]


#r2 = []


'''
for pdb in pdbs:

    #dfols=train_df[train_df['SIZE']==72] # Apenas para 1990
    #dfols=train_df
   
    dfols=train_df[train_df['PDB']!= pdb]
    #dfols=dfols[train_df['PDB']!= 'IT1af7__']
    #dfols=dfols[train_df['PDB']!= 'IT1b72A']
    #dfols=dfols[train_df['PDB']!= 'IT1of9A']
    #dfols=dfols[train_df['PDB']!= 'IT2cr7A']
    #dfols=dfols[train_df['PDB']!= '']
    #dfols= dfols[train_df['PDB']!= pdb]
    
    
    
    #dfols=dfols[train_df['PDB']!= pdb] # Apenas para 1990
    
    
    ols=pd.ols(y=dfols['RMSD']**.3, x=dfols[['SIZE','CONTACT','R_GYRATION','AB_ENERGY','ANGLE','BOND','DIHED','EEL','EELEC','EGB','ESURF','NB','VDWAALS']])

    print '%-12s %14.7f %14.7f ' %(pdb, ols.r2, ols.r2_adj)  
'''



dfols=train_df[train_df['PDB']!= 'IT1kjs_']
#ols=pd.ols(y=dfols['RMSD']**.3, x=dfols[['SIZE','CONTACT','R_GYRATION','R_GYRATION2', 'AB_ENERGY','ANGLE','BOND','DIHED','EEL','EELEC','EGB','ESURF','NB','VDWAALS']])
#ols=pd.ols(y=dfols['RMSD']**.3, x=dfols[keys])
ols=pd.ols(y=dfols['RMSD']**.3, x=dfols[allkeys])
#ols=pd.ols(y=dfols['RMSD']**.3, x=dfols[wkeys])

print (ols)

y = []
yf= []

for i in ols.y:
    y.append(i)

for i in ols.y_predict:
    yf.append(i)


logfile = open('panda_fit.log', 'w')
text = ''

for i in range(0, len(y)):
    text +='%15.7f  %15.7f \n' % ( y[i], yf[i])
    #print text
logfile.write(text)

print ols.beta



