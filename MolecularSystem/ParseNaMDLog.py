'''
ETITLE:      TS           BOND          ANGLE          DIHED          IMPRP               ELECT            VDW       BOUNDARY           MISC        KINETIC               TOTAL           TEMP      POTENTIAL         TOTAL3        TEMPAVG

ENERGY:       0       283.9549       286.8202        56.2013         0.0009            281.1953     20597.2612         0.0000         0.0000       309.6797          21815.1134       295.1475     21505.4338   1129154.1474       295.1475
'''


def SaveNaMDConfigFile(coordinates=None,
                       structure=None,
                       parameters='ff/,par_all27_prot_lipid_na.inp',
                       paratypecharmm='on'):
    """ Function doc """
    pass


def ParseNaMDLog(filein, log=True):
    """ Function doc """
    TS = 0.0
    BOND = 0.0
    ANGLE = 0.0
    DIHED = 0.0
    IMPRP = 0.0
    ELECT = 0.0
    VDW = 0.0
    BOUNDARY = 0.0
    MISC = 0.0
    KINETIC = 0.0
    TOTAL = 0.0
    TEMP = 0.0
    POTENTIAL = 0.0
    TOTAL3 = 0.0
    TEMPAVG = 0.0

    filein = open(filein, 'r')
    for line in filein:
        line2 = line.split()
        if len(line2) > 0:
            if line2[0] == 'ENERGY:':
                # print line
                TS = float(line2[1])
                BOND = float(line2[2])
                ANGLE = float(line2[3])
                DIHED = float(line2[4])
                IMPRP = float(line2[5])
                ELECT = float(line2[6])
                VDW = float(line2[7])
                BOUNDARY = float(line2[8])
                MISC = float(line2[9])
                KINETIC = float(line2[10])
                TOTAL = float(line2[11])
                TEMP = float(line2[12])
                POTENTIAL = float(line2[13])
                TOTAL3 = float(line2[14])
                TEMPAVG = float(line2[15])
    if log == True:
        # print 'TS        = ',TS
        print ('BOND =', BOND,
               'ANGLE =', ANGLE ,
               'DIHED =', DIHED, 
               'IMPRP =', IMPRP, 
               'ELECT =', ELECT, 
               'VDW   =', VDW, 
               'BOUNDARY =', BOUNDARY,
        # print 'MISC      = ',MISC    ,
        # print 'KINETIC   = ',KINETIC ,
        # print 'TEMP      = ',TEMP    ,
              'POTENTIAL =', POTENTIAL)
        # print 'TOTAL     = ',TOTAL
        # print 'TOTAL3    = ',TOTAL3
        # print 'TEMPAVG   = ',TEMPAVG
        # print  BOND + ANGLE + DIHED + IMPRP + ELECT + VDW +  BOUNDARY
        # return POTENTIAL

    enegies = {
        'TS': TS,
        'BOND': BOND,
        'ANGLE': ANGLE,
        'DIHED': DIHED,
        'IMPRP': IMPRP,
        'ELECT': ELECT,
        'VDW': VDW,
        'BOUNDARY': BOUNDARY,
        'MISC': MISC,
        'KINETIC': KINETIC,
        'TOTAL': TOTAL,
        'TEMP': TEMP,
        'POTENTIAL': POTENTIAL,
        'TOTAL3': TOTAL3,
        'TEMPAVG': TEMPAVG,
    }

    #BOND = 0.0
    #ANGLE = 0.0
    #DIHED    = 0.0
    #IMPRP    = 0.0
    #ELECT    = 0.0
    #VDW = VDW * 0.1
    #BOUNDARY = 0.0

    return BOND , ANGLE , DIHED , IMPRP , ELECT , VDW , BOUNDARY


# ParseNaMDLog('HcH_autopsf_MC2.log')
