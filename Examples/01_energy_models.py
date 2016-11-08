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
cmap = CMAP(pdb = os.path.join(PEPDICE_EXAMPLES , 'LABIO_set/1I6C/1I6C_A_AMBER_minimized.pdb'), cutoff = 6.5, log = True)
#-------------------------------------------------------------------------------





# creating a new system 
system = Molecule()
system.name = '1I6C - LABIO dataset' 

# - setup energy model
system.set_energy_model('amber')

# importing coordinates and amber parameters
system.load_PDB_to_system      (filename = os.path.join(PEPDICE_EXAMPLES , 'LABIO_set/1I6C/1I6C_A_AMBER_minimized.pdb'   )   )   
system.import_AMBER_parameters (top      = os.path.join(PEPDICE_EXAMPLES , 'LABIO_set/1I6C/1I6C_A_AMBER.top')   ,   
                                torsions = os.path.join(PEPDICE_PARAMETER, 'amber/AMBER_rotamers.dat') )   


#cmap = CMAP(pdb = os.path.join(PEPDICE_EXAMPLES ,  'Itasser_set/IT1af7__/native_A_AMBER_minimized.pdb'), cutoff = 6.5, log = True)
#system.load_PDB_to_system      (filename = os.path.join(PEPDICE_EXAMPLES , 'Itasser_set/IT1af7__/native_A_AMBER_minimized.pdb') )
#system.import_AMBER_parameters (top      = os.path.join(PEPDICE_EXAMPLES , 'Itasser_set/IT1af7__/native_A_AMBER.top') , 
#                                torsions = os.path.join(PEPDICE_PARAMETER, 'amber/AMBER_rotamers.dat') )   
#



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

