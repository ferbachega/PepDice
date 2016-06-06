import subprocess


#---------------------------------------------------------------------------------------------------#
import os                                                                                           #
from pprint      import pprint                                                                      #
from Molecule    import Molecule                                                                    #
from Geometry    import *                                                                           #
from MonteCarlo  import monte_carlo,  monte_carlo_dic, MC_replica_exchange, run_MC_replica_exchange #
from XYZFiles    import save_XYZ_to_file                                                            #
from CRDFiles    import load_CRD_from_file                                                          #
from AATorsions  import ROTAMER_LIST                                                                #
                                                                                                    #
from RMSD import compute_RMSD                                                                       #
                                                                                                    #
from Test        import *                                                                           #
                                                                                                    #
from GeometryOptimization import minimize                                                           #
from Energy import save_PDB_to_file                                                                 #
                                                                                                    #
#---------------------------------------------------------------------------------------------------#
from SideChainRefine import optimize_side_chain                                                     #
#---------------------------------------------------------------------------------------------------#
                                                                                                    #

PEPDICE = os.environ.get('PEPDICE')
PEPDICE_EXAMPLES = os.path.join(PEPDICE, 'Examples')
PEPDICE_PARAMETER= os.path.join(PEPDICE, 'Parameters')







def  build_tleap( pdbin , ff):
    """ Function doc """
    
    text  = 'source leaprc.' + ff + ' \n'
    text += 'foo = loadpdb ' + pdbin       + ' \n'
    text += 'saveamberparm foo '+ pdbin[:-4] +ff+'.top '+  pdbin[:-4]+ ff+'.crd \n'
    text += 'savepdb foo ' + pdbin[:-4]+ ff+'.pdb \n'
    text += 'quit'
    leaprc = open('leaprc', 'w')
    leaprc.write(text)
    leaprc.close()
    subprocess.call(['tleap', '-f', 'leaprc'])
    return pdbin[:-4] + ff








pdbs = ['2n1p',
        '2n2r',
        '2n2s',
        '2n52',
        '2n5j',
        '2n5u',
        '2n7f',
        '2n7i',
        '2n9c',
        '2nat',
        '2nav',
        '2nbd',
        '2nbe']




for pdb in pdbs:
    
    text  = ''
    model = 1 
    pdbfile = open(pdb + '.pdb', 'r')
    for line in pdbfile:
        
        line2 = line.split()
        if line2[0] ==  'ENDMDL':
            
            pdbout = open(pdb +'_'+str(model)+'_'+'.pdb', 'w')
            pdbout.write(text)
            pdbout.close()

            pdbout = build_tleap( pdb +'_'+str(model)+'_'+'.pdb', 'ff03ua')


            system = Molecule() 
            system.load_PDB_to_system         (filename = pdbout +'.pdb')   
            system.import_AMBER_parameters    (top      = pdbout +'.top', 
                                               torsions = os.path.join(PEPDICE_PARAMETER, 'amber/AMBER_rotamers.dat')) 

            
            print system.energy()
            minimize(molecule = system,
                     imin  = 1          ,
                     maxcyc= 5000       ,
                     ncyc  = 2000       ,
                     cut   = 12         ,
                     rgbmax= 999        ,
                     igb   = 1          ,
                     ntb   = 0          ,
                     ntpr  = 100        ,
                     ntr   = 0          )



            print system.energy()
            save_PDB_to_file(system,  pdbout + '_minimized.pdb')






            model += 1 
            text = ''
        


        else:
            if line2[0] ==  'ATOM' and line2[-1] != 'H':
                text += line
            else:
                pass
                


