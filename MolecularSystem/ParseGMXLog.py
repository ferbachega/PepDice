def ParseGMXLog(filein, log=True):
  
    keys = {
            'Angle'          : None,
            'Proper Dih.'    : None,
            'Improper Dih.'  : None,
            'LJ-14'          : None,
            'Coulomb-14'     : None,
            'Potential'      : None,
            }
    

    
    
    filein = open(filein, 'r')
    arq2   = filein.readlines()
    
    for line in arq2:
        line2 = line.split()
     
        if len(line2) == 8:
            keys['Angle'        ] = line2[1]
            keys['Proper Dih.'  ] = line2[2]
            keys['Improper Dih.'] = line2[3]
            keys['LJ-14'        ] = line2[4]
            keys['Coulomb-14'   ] = line2[5]
            keys['Potential'    ] = line2[7]
    if log:
        from pprint import pprint
        pprint(keys)
    return keys
