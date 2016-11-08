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



text = '%-20s %-10s %6s  %6s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s' %('DECOY',
                                                                                        'PDB',
                                                                                        'RMSD',
                                                                                        'SIZE',
                                                                                        'AB_ENERGY',
                                                                                        'CONTACT',
                                                                                        'ANGLE',
                                                                                        'BOND',
                                                                                        'DIHED',
                                                                                        'EEL',
                                                                                        'EELEC',
                                                                                        'EGB',
                                                                                        'ESURF',
                                                                                        'NB',
                                                                                        'VDWAALS')  

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
    system.set_energy_model('amber')
    system.import_AMBER_parameters (top      = 'native_A_AMBER.top'                ,   
                                    torsions = os.path.join(PEPDICE_PARAMETER, 'amber/AMBER_rotamers.dat') )       
    
    cmap = CMAP(pdb = 'native_A_AMBER_minimized.pdb', cutoff = 6.5, log = False)
    
    
    
    
    for pdb in RMSD_list:
        if pdb == 'NAME':
            pass
        
        else:
            filename =  pdb.replace('.', '_A_AMBER_minimized.')
            system.load_PDB_to_system      (filename = filename)   
            system.set_energy_model('FULL')

            #pprint(RMSD_list)
            system.import_CMAP(cmap = cmap)
            
            try:
                energies = system.energy(return_list = True)
                
                
                
                text = '%-20s %-10s %6s %6d %15.7f %15.7f %15.7f %15.7f %15.7f %15.7f %15.7f %15.7f %15.7f %15.7f %15.7f' %( pdb                   ,
                                                                                                                             folder                ,
                                                                                                                             RMSD_list[pdb]        ,
                                                                                                                             len(system.residues)  , 
                                                                                                                              energies['AB_ENERGY'],
                                                                                                                              energies['CONTACT'  ],
                                                                                                                              energies['ANGLE'    ],
                                                                                                                              energies['BOND'     ],
                                                                                                                              energies['DIHED'    ],
                                                                                                                              energies['EEL'      ],
                                                                                                                              energies['EELEC'    ],
                                                                                                                              energies['EGB'      ],
                                                                                                                              energies['ESURF'    ],
                                                                                                                              energies['NB'       ],
                                                                                                                              energies['VDWAALS'  ],
                                                                                                                              )
                
                if 'nan' in text:
                    pass
                else:
                    textlines.append(text)
                    logfile.writelines(text)
                    logfile.write('\n')    
                    print text                                                                                           
            except:
                pass
            
 
