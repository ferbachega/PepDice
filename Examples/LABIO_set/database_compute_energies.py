# ------------------------------------------------------------------------------
import os
import subprocess
import tempfile

from Energy import save_PDB_to_file
from GeometryOptimization import minimize
from Molecule import Molecule
from CMAP import CMAP
# ------------------------------------------------------------------------------
from pprint import pprint

PEPDICE = os.environ.get('PEPDICE')
PEPDICE_EXAMPLES = os.path.join(PEPDICE, 'Examples')
PEPDICE_PARAMETER = os.path.join(PEPDICE, 'Parameters')

PEPDICE_PDBS = os.path.join(PEPDICE, 'PDBs')
PEPDICE_OUTPUTS = os.path.join(PEPDICE, 'outputs')


def  import_rmsd_from_file ( filein = None):
    """ Function doc """
    text =  open(filein , 'r')
    
    RMSD_list = {}
    
    for line in text:
        line2 = line.split()
        
        if line2 == 'NAME':
            pass
        
        else:
            if len(line2) == 2:
               #print line2[1]
               RMSD_list[line2[0]] =line2[1]

    return RMSD_list


cwd = os.getcwd()

folders = [
        '1BDD',
        '1CSK',
        '1DV0',
        '1E0L',
        '1E0N',
        '1EI0',
        '1ERY',
        '1FME',
        '1GAB',
        '1I6C',
        '1L2Y',
        '1PRB',
        '1PSV',
        '1RES',
        '1RIJ',
        '1VII',
        '1WY3',
        '1YRF',
        '2F21',
        '2F4K',
        '2HBA',
        '2HEP',
        '2JOF',
        '2WXC',
        ]



text = '%-20s %-10s %6s  %6s %15s ' %('DECOY',
                                      'PDB',
                                      'RMSD',
                                      'SIZE',
                                      'energy',)  


textlines = []
textlines.append(text)
logfile =  open('logfile.txt', 'w')
logfile.writelines(text)
logfile.write('\n')

print text


for folder in folders:
    path = os.path.join(cwd, folder)
    
    os.chdir(os.path.join(cwd, folder))
    
    RMSD_list = import_rmsd_from_file ( filein = 'list.txt')
    system = Molecule()
    system.set_energy_model('LSF')
    system.import_AMBER_parameters (top      = folder+'_A_AMBER.top'                ,   
                                    torsions = os.path.join(PEPDICE_PARAMETER, 'amber/AMBER_rotamers.dat') )       
    
    
    
    
    
    
    for pdb in RMSD_list:
        if pdb == 'NAME':
            pass
        
        else:
            filename =  pdb.replace('.', '_A_AMBER_minimized.')
            system.load_PDB_to_system      (filename = filename)   
            system.set_energy_model('LSF')

            #pprint(RMSD_list)
            

            
            try:    
                cutoff   = 6.0
                contacts = []
                
                
                #saltcon = 0.00
                #EGBs = []
                #
                #for i in range(0, 6):
                #    
                #    system.energy_model_parameters['saltcon'] = saltcon
                #    energies = system.energy(return_list = True)
                #    saltcon += 0.1
                #    EGBs.append(energies['EGB'])
                
                #energies = system.energy(return_list = True)
                
                energy = system.energy()
                
                #for i in range (0, 6):  
                #    cmap = CMAP(pdb = 'native_A_AMBER_minimized.pdb', cutoff = cutoff, log = False)
                #    system.import_CMAP(cmap = cmap)
                #    contact =  system.compute_CONTACT_energy()
                #    contacts.append(contact)
                #    cutoff += 0.5
                
                #print folder
                #print contacts
                

                text = '%-20s %-10s %6s %6d %15.7f ' %( pdb                   ,
                                                        folder                ,
                                                        RMSD_list[pdb]        ,
                                                        len(system.residues)  ,
                                                        energy                ,
                                                        )

    
                textlines.append(text)
                
                logfile.writelines(text)
                logfile.write('\n')    
                
                print text                                                                                           
            except:
                pass

 
