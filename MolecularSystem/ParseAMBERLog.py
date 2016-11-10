'''
 NSTEP  =        0        TIME(PS) =       0.000    TEMP(K)    =   299.72       PRESS =     0.0
 Etot   =    -53913.3910  EKtot    =     23324.7176  EPtot      =    -77238.1086
 BOND   =       380.2130  ANGLE    =       934.4207  DIHED      =      2366.0987
 1-4 NB =      1231.0477  1-4 EEL  =     12302.5221  VDWAALS    =      8432.2318
 EELEC  =    -87295.5824  EGB      =    -16446.7036  RESTRAINT  =         0.0000
 ESURF  =       857.6436'''


def SaveNaMDConfigFile(coordinates=None,
                       structure=None,
                       parameters='ff/,par_all27_prot_lipid_na.inp',
                       paratypecharmm='on'):
    """ Function doc """
    pass


def ParseAMBERLog(filein, log=True):
    """ Function doc """
    keys = { #"NSTEP": None, 
          #"TIME(PS)": None,
           #"TEMP(K)": None,
             #"PRESS": None,
              #'Etot': None,
             #'EKtot': None,
             #'EPtot': None,
              'BOND': None,
             'ANGLE': None,
             'DIHED': None,
                'NB': None,
               'EEL': None,
           'VDWAALS': None,
             'EELEC': None,
               'EGB': None,
         #'RESTRAINT': None,
             'ESURF': None, 
           'CONTACT': 0.0 ,
         'AB_ENERGY': 0.0 ,    
             }  #'ESURF=' ]    
    
    
    NSTEP     = 0
    Etot      = 0
    BOND      = 0
    NB14      = 0
    EELEC     = 0
    ESURF     = 0
    TIME      = 0
    EKtot     = 0
    ANGLE     = 0
    EEL14     = 0
    EGB       = 0
    TEMP      = 0
    EPtot     = 0
    DIHED     = 0
    VDWAALS   = 0
    RESTRAINT = 0
    PRESS     = 0
    
    
    filein = open(filein, 'r')
    arq2   = filein.readlines()
    
    for line in arq2:
        line = line.replace('=','')
        line2 = line.split()
     
        for string in line2:
            if string in keys:
                try:
                    keys[string] = float(line2[line2.index(string)+1])
                except:
                    pass
    return keys
    #from pprint import pprint
    #pprint (keys)
                #print string, line2[line2.index(string)+1]
                
                
                

#ParseAMBERLog('/home/farminf/Programas/PepDice/Examples/1GAB/energy.log')
