#----------------------------------------------------#
import os                                            #
from Molecule    import Molecule                     #
#----------------------------------------------------#

#-------------------------------------------------------------------------------
PEPDICE = os.environ.get('PEPDICE')
PEPDICE_EXAMPLES = os.path.join(PEPDICE, 'Examples')
PEPDICE_PARAMETER= os.path.join(PEPDICE, 'Parameters')
#-------------------------------------------------------------------------------


# creating a new system 
system = Molecule()
system.name = '1I6C - LABIO dataset' 

# - energy model AMBER
system.set_energy_model('amber')

# importing coordinates and amber parameters
system.load_PDB_to_system      (filename = os.path.join(PEPDICE_EXAMPLES , 'LABIO_set/1I6C/1I6C_A_AMBER_minimized.pdb'   )   )   
system.import_AMBER_parameters (top      = os.path.join(PEPDICE_EXAMPLES , 'LABIO_set/1I6C/1I6C_A_AMBER.top')   ,   
                                torsions = os.path.join(PEPDICE_PARAMETER, 'amber/AMBER_rotamers.dat') )   

system.Status()




