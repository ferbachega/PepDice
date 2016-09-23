#-------------------------------------------------------------------------------
import os                                             
from pprint      import pprint                        
from Molecule    import Molecule                      
from Geometry    import *                             
from MonteCarlo  import monte_carlo, insert_fragment  
from XYZFiles    import save_XYZ_to_file              
from CRDFiles    import load_CRD_from_file            
from AATorsions  import ROTAMER_LIST                  
                                                      
from RMSD import compute_RMSD                         
                                                                                                           
from GeometryOptimization import minimize             
                                                      
from random import randint                            
                                                      
from Energy import save_PDB_to_file                   
#from Amber12ToAmber11 import amber12_to_amber11_topology_converter

#-------------------------------------------------------------------------------
from Fragments import import_fragments_from_pdb
import subprocess
#-------------------------------------------------------------------------------

PEPDICE = os.environ.get('PEPDICE')
PEPDICE_EXAMPLES = os.path.join(PEPDICE, 'Examples')
PEPDICE_PARAMETER= os.path.join(PEPDICE, 'Parameters')

PEPDICE_PDBS     = os.path.join(PEPDICE, 'PDBs')
PEPDICE_OUTPUTS  = os.path.join(PEPDICE, 'outputs')


os.path.join(PEPDICE_PDBS,)
#-------------------------------------------------------------------------------


def pdb_hydrogen_remove (pdbin = None, pdbout = False):
    """ Function doc """
    print pdbin
    print pdbout
    

def pdb_extract_chain (pdbin = None, pdbout = None, chain = 'A', model = 1, remove_hydrogens = True):
    """ Function doc """
    
    pdbtext = open(pdbin, 'r')
    text    = ''
    for line in pdbtext:
        line2 = line.split()
        #print line
        if len(line2) > 0:
            if line2[0] == 'ATOM':
            
                if chain in line2:
                    
                    # ignora as linhas com H no final
                    if line2[-1] == 'H':
                        pass
                    
                    
                    else:
                        try:
                            text += line +'\n' 
                        except:
                            print line
                else:
                    pass
            
            if line2[0] == 'MODEL':
                if int(line2[1]) == 1:
                    pass
                else:
                    break
                
    pdbout = open(pdbout, 'w')
    pdbout.write(text)
    pdbout.close()


def get_sequence_from_pdb (pdbin = None, seq_out = None):
    """ Function doc """
    pdb     = pdbin
    
    pdbtext =  open(pdb, 'r')
    sequence      = []
    sequence_code = ''
    previous_index = None

    aa_dic = { 
         'ALA':'A',
         'ARG':'R',
         'ASN':'N',
         'ASP':'D',
         'CYS':'C',
         'GLU':'E',
         'GLN':'Q',
         'GLY':'G',
         'HIS':'H',
         'ILE':'I',
         'LEU':'L',
         'LYS':'K',
         'MET':'M',
         'PHE':'F',
         'PRO':'P',
         'SER':'S',
         'THR':'T',
         'TRP':'W',
         'TYR':'Y',
         'VAL':'V'
         }
    
    
    for line in pdbtext:
        line2 = line.split()
        #print line
        if len(line2) > 0:
            if line2[0] == 'ATOM':
                #print line
                
                try:
                    index = int(line2[5])
                    
                    if index == previous_index:
                        pass
                    else:
                        sequence.append(line2[3])
                        previous_index = index
                        sequence_code += aa_dic[line2[3]]
                except:
                    print line
                #print index
                #sequence[index] = line2[3]
                
            if line2[0] == 'TER' or line2[0] == 'END' :
                break
    
    if seq_out :
        pdbcode = pdbin.split('.')
        pdbcode = pdbcode[0]
        pdbcode = pdbcode +'  ' + sequence_code
        fileout = open(seq_out, 'w')
        fileout.write(pdbcode)
    
    return sequence_code


def phenix_geometry_minimization (pdbin = None, geo = False, geofile = None ):
    """ Function doc """
    
    if geo:
        subprocess.call(['phenix.geometry_minimization', pdbin, geofile])
    else:
        subprocess.call(['phenix.geometry_minimization', pdbin])
    
    return pdbin+'_minimized.pdb'


def build_AMBER_system_from_PDB (pdbin = None    ,
                              basename = None    ,
                           force_field = 'ff03ua',
                             overwrite = True 
                                 ):
    """ 
    Function doc 
    source leaprc.ff03ua
    foo = sequence { ACE ALA NME }
    saveamberparm foo foo.top foo.crd
    """
    
    if overwrite:
        #text  = 'source leaprc.' + force_field + ' \n'
        text  = 'addpath '+PEPDICE+'/Parameters/amber/labio.amber \n'
        text  += 'source leaprc.' + force_field + ' \n'

        text += 'foo = loadpdb '+ pdbin
        
        text += '\n'
        text += 'saveamberparm foo '+ basename +'.top '+ basename +'.crd \n'
        text += 'savepdb foo ' + basename +'.pdb \n'
        text += 'quit'
        leaprc = open('leaprc', 'w')
        leaprc.write(text)
        leaprc.close()

    subprocess.call(['tleap', '-f', 'leaprc'])
    #amber12_to_amber11_topology_converter (basename +'.top', basename +'.top')


def amber12_to_amber11_topology_converter (filein, fileout):
	filein = open(filein, 'r')
	text   = []
	print_line = True

	for line in filein:
		line2 = line.split()
		try:
			if line2[0] == '%FLAG':
				if   line2[1] == 'ATOMIC_NUMBER':
					print 'excluding flag:', line
					print_line = False

				elif   line2[1] == 'SCEE_SCALE_FACTOR':
					print 'excluding flag:', line
					print_line = False

				elif   line2[1] == "SCNB_SCALE_FACTOR":
					print 'excluding flag:', line
					print_line = False			

				elif   line2[1] == 'IPOL':
					print 'excluding flag:', line
					print_line = False
		
				else:
					print_line = True	
					#print print_line
		except:
			a= None
		if print_line == True:
			text.append(line)

	fileout = open(fileout, 'w')
	fileout.writelines(text)
	fileout.close()


def amber_topology_angle_force_change (filein =  None, fileout = None, force = 'E+04'):
    """ Function doc """
    filein =  open(filein, 'r')
    filein2 = filein.readlines()

    text = []
    
    for line in filein2:
        if '%FLAG ANGLE_FORCE_CONSTANT                                                      ' in line:
            inicio = filein2.index(line)

        if '%FLAG ANGLE_EQUIL_VALUE                                                         ' in line:
            final = filein2.index(line)

    #print inicio, final
    for line in range(inicio, final):
        print filein2[line]
        
        if 'E+01' in filein2[line]:
            filein2[line] = filein2[line].replace('E+01', force)
            print filein2[line]


    fileout =  open(fileout, 'w')
    fileout.writelines(filein2)
    fileout.close()
    

def wget_pdb (pdbcode = None, filesInFolder = None):
    """ Function doc """
    
    if filesInFolder == None:
        filesInFolder = os.listdir('.')
    
    if pdbcode+'.pdb' in filesInFolder:
        pass
    else:
        os.system('wget https://files.rcsb.org/download/'+pdbcode+'.pdb') 


PDBcodes =  {
            #'1GAB': ['alpha'     ],
            #'1BX4': ['alpha_beta'],
             #'1UAO': ['beta'      ],
             #'1LE1': ['beta'      ],
             #'1CSK': ['beta'      ],
             #'2GB1': ['beta'      ],
             #'1E0Q': ['beta'      ],
             #'1L2Y': ['alpha'     ],
             #'2JOF': ['alpha'     ],
             #'1RIJ': ['alpha'     ],
             #'1E0N': ['beta'      ],
             #'1E0L': ['beta'      ],
             #'1I6C': ['beta'      ],
             #'1FME': ['alpha_beta'],
             #'1PSV': ['alpha_beta'],
             #'1UBQ': ['alpha_beta'],
             #'1WY3': ['alpha'     ],
             #'1YRF': ['alpha'     ],
             #'2F4K': ['alpha'     ],
             #'1VII': ['alpha'     ],
             #'1EI0': ['alpha'     ],
             #'1ERY': ['alpha'     ],
             #'2HBA': ['alpha_beta'],
             #'2HEP': ['alpha'     ],
             #'1RES': ['alpha'     ],
             #'1BDD': ['alpha'     ],
             #'1E0G': ['alpha_beta'],
             #'1BDD': ['alpha'     ],
             '1DV0': ['alpha'     ],
             '1PRB': ['alpha'     ],
             '2WXC': ['alpha'     ],
            }




folder = os.getcwd()
Types = ['alpha', 'beta', 'alpha_beta']

for Type in Types:
    if not os.path.exists (Type):
        os.mkdir (Type)

logfile = open('logfile', 'w')
logtext = []
for code in PDBcodes:
    print code
    print PDBcodes[code][0]
    #---------------------------------------------------------------------------
    path = os.path.join(folder, PDBcodes[code][0])
    os.chdir(path)
    
    
    if not os.path.exists (code):
        os.mkdir (code)
    
    os.chdir(os.path.join(path, code))
    #---------------------------------------------------------------------------
    
    
    
    try:
        #---------------------------------------------------------------------------
        wget_pdb (pdbcode = code, filesInFolder = None)
        
        pdb_extract_chain (pdbin = code+'.pdb', 
                          pdbout = code+'_A.pdb', 
                           chain = 'A')
        

        get_sequence_from_pdb      (pdbin = code+'_A.pdb', seq_out = code+'_A.seq')    
        #---------------------------------------------------------------------------
        logtext.append(code+' wget_pdb......Ok\n')
    except:
        print 'failed wget pdb'
        logtext.append(code+' wget_pdb......failed\n')

    
    
    try:
        #---------------------------------------------------------------------------
        build_AMBER_system_from_PDB(pdbin = code+'_A.pdb'    ,
                                 basename = code+'_A_AMBER'    ,
                              force_field = 'ff03ua.labio'         ,        
                                overwrite = True 
                                    )
        #----------------------------------------------------------------------------
        
        
        ##----------------------------------------------------------------------------
        #amber_topology_angle_force_change(filein = code+'_A_AMBER.top', 
        #                                 fileout = code+'_A_AMBER_angleMod.top', 
        #                                   force = 'E+04')
        ##----------------------------------------------------------------------------
        
        
        #----------------------------------------------------------------------------
        system = Molecule() 
        system.load_PDB_to_system      (filename = code+'_A_AMBER.pdb')   
        system.import_AMBER_parameters (top      = code+'_A_AMBER.top',   
                                        torsions = os.path.join(PEPDICE_PARAMETER, 'amber/AMBER_rotamers.dat') )
        save_PDB_to_file(system, code+'_A.pdb')
        logtext.append(code+' tleap......Ok\n')

    except:
        print 'failed wget pdb'
        logtext.append(code+' tleap......failed\n')

    #'''
    try:
        minimize(molecule = system,
                 imin  = 1        ,
                 maxcyc= 1000     ,
                 ncyc  = 100      ,
                 cut   = 10       ,
                 rgbmax= 999      ,
                 igb   = 1        ,
                 ntb   = 0        ,
                 ntpr  = 100      ,
                 ntr   = 0        )
                 #restraintmask = ':1-50 & @CA,N,C,O=', 
                 #restraint_wt  =  50.0 
        save_PDB_to_file(system, code+'_A_AMBER_minimized.pdb')
        logtext.append(code+' amber opt......ok\n')
    
        #----------------------------------------------------------------------------
    except:
        print 'failed:', code, 'opt amber'
        logtext.append(code+' amber opt......failed\n')
    #'''


logfile.writelines(logtext)



'''
get_sequence_from_pdb (pdbin = '1GAB_noH.pdb', seq_out = '1GBA.seq')    

build_AMBER_system_from_PDB(pdbin = '1GAB_noH.pdb'    ,
                         basename = '1GAB_noH'        ,
                      force_field = 'ff03ua',
                        overwrite = True 
                            )
'''