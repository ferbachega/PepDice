from ParseNaMDLog import ParseNaMDLog
from ParseAMBERLog import ParseAMBERLog
from CRDFiles      import save_CRD_to_file, load_CRD_from_file
import os
import sys


def minimize(molecule  = None                        ,
                   log = False                       , 
                 imin  = 1                           ,
                 maxcyc= 1000                        ,
                 ncyc= 100                           ,
                 cut=10                              ,
                 rgbmax=999                          ,
                 igb=1                               ,
                 ntb=0                               ,
                 ntpr=100                            ,
                 ntr = 0                             ,
                 restraintmask = ':1-50 & @CA,N,C,O=', 
                 restraint_wt  =  50.0               , 
              ):

    
    if  molecule.ff_type == 'amber':
        
        write_AMBER_minimization_input_file(molecule      =molecule      , 
                                            filename      ='minimize.in' ,
                                            imin          =imin          ,
                                            maxcyc        =maxcyc        ,
                                            ncyc          =ncyc          ,
                                            cut           =cut           ,
                                            rgbmax        =rgbmax        ,
                                            igb           =igb           ,
                                            ntb           =ntb           ,
                                            ntpr          =ntpr          ,
                                            ntr           =ntr           ,
                                            restraintmask =restraintmask , 
                                            restraint_wt  =restraint_wt )
        
        
        save_CRD_to_file (molecule = molecule, filename='tmp.crd')

        os.system('sander -O -i minimize.in -c tmp.crd -o minimized.log -p ' + molecule.top + ' -r minimized.rst -ref tmp.crd')
        load_CRD_from_file(molecule = molecule, filename = 'minimized.rst')











def write_AMBER_minimization_input_file (molecule =None                      , 
                                         filename = 'minimize.in'            ,
                                         imin  = 1                           ,
                                         maxcyc= 1000                        ,
                                         ncyc= 100                           ,
                                         cut=10                              ,
                                         rgbmax=999                          ,
                                         igb=1                               ,
                                         ntb=0                               ,
                                         ntpr=100                            ,
                                         ntr = 0                             ,
                                         restraintmask = ':1-50 & @CA,N,C,O=', 
                                         restraint_wt  =  50.0               ):
                                             
                                             
    text = """Stage 1 - minimisation of TC5b
 &cntrl
  imin=%i, maxcyc=%i, ncyc=%i,
  cut=%i, rgbmax=%i,igb=%i, ntb=%i,
  ntpr=%i, 
  ntr=%i,
  restraintmask = '%s',
  restraint_wt=%3.1f,
 /
 """ %(imin,maxcyc,ncyc, cut,rgbmax,igb,ntb,ntpr, ntr,restraintmask, restraint_wt)
 
    output_file = open(filename, "w")
    output_file.write(text)
    output_file.close()
