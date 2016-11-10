#----------------------------------------------------#
import os                                            #
from pprint      import pprint                       #
from Molecule    import Molecule                     #
#----------------------------------------------------#

#                         environment variables 
#-------------------------------------------------------------------------------
PEPDICE = os.environ.get('PEPDICE')
PEPDICE_EXAMPLES = os.path.join(PEPDICE, 'Examples')
PEPDICE_PARAMETER= os.path.join(PEPDICE, 'Parameters')
#-------------------------------------------------------------------------------



# building a contact map - required for Contact model calculations
#-------------------------------------------------------------------------------
from CMAP import CMAP
cmap = CMAP(pdb = os.path.join(PEPDICE_EXAMPLES , 'LABIO_set/1GAB/1GAB_A_AMBER_minimized.pdb'), cutoff = 6.5, log = True)
#-------------------------------------------------------------------------------





# creating a new system 
system = Molecule()
system.name = '1GAB - LABIO dataset' 

# - setup energy model
system.set_energy_model('amber')

# importing coordinates and amber parameters
system.load_PDB_to_system      (filename = os.path.join(PEPDICE_EXAMPLES , 'LABIO_set/1GAB/1GAB_A_AMBER_minimized.pdb'   )   )   
system.import_AMBER_parameters (top      = os.path.join(PEPDICE_EXAMPLES , 'LABIO_set/1GAB/1GAB_A_AMBER.top')   ,   
                                torsions = os.path.join(PEPDICE_PARAMETER, 'amber/AMBER_rotamers.dat') )   



print system.compute_R_gy_Calpha()
#print system.compute_SS_energy(log = True)
#cmap = CMAP(pdb = os.path.join(PEPDICE_EXAMPLES ,  'Itasser_set/IT1af7__/native_A_AMBER_minimized.pdb'), cutoff = 6.5, log = True)
#system.load_PDB_to_system      (filename = os.path.join(PEPDICE_EXAMPLES , 'Itasser_set/IT1af7__/native_A_AMBER_minimized.pdb') )
#system.import_AMBER_parameters (top      = os.path.join(PEPDICE_EXAMPLES , 'Itasser_set/IT1af7__/native_A_AMBER.top') , 
#                                torsions = os.path.join(PEPDICE_PARAMETER, 'amber/AMBER_rotamers.dat') )   
#


system.import_SS_restraints_from_string (      ss = 'CCHHHHHHHHHHHHHHHHHHCCCCCHHHHHHHHHCCCHHHHHHHHHHHHHHCC', 
                                            w_ss  = '00000123345553221111000000001234432100012342335556110', log= True)
print system.compute_SS_energy(log = True)

#restraint = {
#                    'Kd'            : 40.0,
#                    'resi_i'        : 1   , #starts at ZERO
#                    'resi_j'        : 5   , 
#                    'atom_name_i'   : 'CA',
#                    'atom_name_j'   : 'CA',
#                    'distance'      : 10.0,
#                    }
#
#system.hamonical_potential_restraint_list.append(restraint)

print system.compute_harmonical_restraint_energies(log = True)

energy_models = ['LABIO']#, 'FULL', 'iLABIO']


for model in energy_models:
    system.set_energy_model(model)
    #system.energy_model_parameters['saltcon'] = 1.0
    #system.energy_model_parameters['surften'] = 1.0
    #system.energy_model_parameters['igb'] = 1
    #if model == 'Contact':
    system.import_CMAP(cmap = cmap)
    
    system.Status()
    system.energy( log =True)
    #print 'energy:', model ,system.energy()
    print '\n\n\n'

