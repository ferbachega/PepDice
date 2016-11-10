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



from Geometry    import  distance_ab, computePhiPsi



residue_type_dic = {'ALA': [],
                    'ARG': [],
                    'ASN': [],
                    'ASP': [],
                    'CYS': [],
                    'GLU': [],
                    'GLN': [],
                    'GLY': [],
                    'HIS': [],
                    'ILE': [],
                    'LEU': [],
                    'LYS': [],
                    'MET': [],
                    'PHE': [],
                    'PRO': [],
                    'SER': [],
                    'THR': [],
                    'TRP': [],
                    'TYR': [],
                    'VAL': [],
                    }


for folder in folders:
    path = os.path.join(cwd, folder)
    
    os.chdir(os.path.join(cwd, folder))
    system = Molecule()
    system.set_energy_model('FULL')
    
    
    system.import_AMBER_parameters (top      = folder+'_A_AMBER.top'                ,   
                                    torsions = os.path.join(PEPDICE_PARAMETER, 'amber/AMBER_rotamers.dat') )       
    
    system.load_PDB_to_system      (filename = 'native_A_AMBER_minimized.pdb')
    
    R_abs = []
    
    
    
    
    for i in range(0, len(system.residues)):
        for j in range(i+1 , len(system.residues)):
            
            name_i = system.residues[i].name
            for atom in system.residues[i].atoms:
                if atom.name == 'CA':
                    atom_i    = atom  
                    atom_i.AB          = 'A'
                    
                    

            name_j = system.residues[j].name
            for atom in system.residues[j].atoms:
                if atom.name == 'CA':
                    atom_j = atom  
                    atom_j.AB = 'A'            
            
            R_ab = distance_ab (atom_i, atom_j)
            #print R_ab
            R_abs.append(R_ab)
            
            if R_ab:
            
                #if type(R_ab) == float:
            
                    residue_type_dic[system.residues[i].name].append(R_ab)
                    #print system.residues[i].name, R_ab , system.residues[j].name
                #else:
                #    print 'R_ab is not a float'
            
        
from pprint import pprint
pprint(R_abs)

#for key in residue_type_dic:
#    print key
#    pprint (residue_type_dic[key])




#for i in range(0, len(residue_type_dic['ALA'])):
#    print(residue_type_dic['ALA'][i],
#          residue_type_dic['ARG'][i],
#          residue_type_dic['ASN'][i],
#          residue_type_dic['ASP'][i],
#          residue_type_dic['CYS'][i],
#          residue_type_dic['GLU'][i],
#          residue_type_dic['GLN'][i],
#          residue_type_dic['GLY'][i],
#          residue_type_dic['HIS'][i],
#          residue_type_dic['ILE'][i],
#          residue_type_dic['LEU'][i],
#          residue_type_dic['LYS'][i],
#          residue_type_dic['MET'][i],
#          residue_type_dic['PHE'][i],
#          residue_type_dic['PRO'][i],
#          residue_type_dic['SER'][i],
#          residue_type_dic['THR'][i],
#          residue_type_dic['TRP'][i],
#          residue_type_dic['TYR'][i],
#          residue_type_dic['VAL'][i])
#
