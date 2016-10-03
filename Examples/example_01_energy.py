#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  example01.py
#  
#  Copyright 2016 Fernando Bachega <fernando@Fenrir>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  


#----------------------------------------------------#
import os                                            #
from pprint      import pprint                       #
from Molecule    import Molecule                     #
from Geometry    import *                            #
from MonteCarlo  import monte_carlo                  #
from XYZFiles    import save_XYZ_to_file             #
from CRDFiles    import load_CRD_from_file           #
from AATorsions  import ROTAMER_LIST                 #
                                                     #
from RMSD import compute_RMSD                        #
                                                     #
from Test        import *                            #
                                                     #
from GeometryOptimization import minimize            #
                                                     #
from random import randint                           #
                                                     #
from Energy import save_PDB_to_file                  #
                                                     #
#----------------------------------------------------#
from SideChainRefine import optimize_side_chain
#----------------------------------------------------#
from Fragments import build_fragment_library_from_pdbs














#eita

#-------------------------------------------------------------------------------
PEPDICE = os.environ.get('PEPDICE')
PEPDICE_EXAMPLES = os.path.join(PEPDICE, 'Examples')
PEPDICE_PARAMETER= os.path.join(PEPDICE, 'Parameters')
#-------------------------------------------------------------------------------











#-------------------------------------------------------------------------------
system = Molecule() 
system.load_PDB_to_system      (filename = os.path.join(PEPDICE_EXAMPLES , 'PDBs/alpha/1GAB/1GAB_A_AMBER_minimized.pdb'   )   )   
system.import_AMBER_parameters(top       = os.path.join(PEPDICE_EXAMPLES , 'PDBs/alpha/1GAB/1GAB_A_AMBER.top')   ,   
                                torsions = os.path.join(PEPDICE_PARAMETER, 'amber/AMBER_rotamers.dat') )   
#-------------------------------------------------------------------------------
system.AB = 1
print system.energy(log       = True,
                    AB_energy = True)










#-------------------------------------------------------------------------------
system = Molecule() 
system.load_PDB_to_system      (filename = os.path.join(PEPDICE_EXAMPLES , 'PDBs/alpha/1GAB/1GAB_A_AMBER_minimized.pdb'   )   )   

system.ff_type = 'Calpha_model'
system.import_Calpha_model_parameters(top       = os.path.join(PEPDICE_EXAMPLES , 'PDBs/alpha/1GAB/1GAB_A_AMBER.top')   ,   
                                      torsions  = os.path.join(PEPDICE_PARAMETER, 'amber/AMBER_rotamers.dat') )  
TRAJECTORY = 'Calpha_test'                                      
                                      

#-------------------------------------------------------------------------------
system.AB = 1
print system.energy(log       = True,
                    AB_energy = True)

















'''
#--------------------------------------------------------------------------------------------------------------------------------------------

files =  os.listdir('/home/fernando/programs/pepdice/Examples/teste/RO1aiu/')

rmsd  =  {}

text = open('/home/fernando/programs/pepdice/Examples/teste/RO1aiu/' + 'list.txt', 'r')
for line in text:
    line2 = line.split()
    if len(line2) == 2:
        if line2[0] == 'NAME':
            pass
        else:
            rmsd[line2[0]] = line2[1]



for _file in files:
    file2 = _file.split('.')
    
    #print file2
    
    if file2[-1] == 'pdb':
        
        system = Molecule() 
        system.load_PDB_to_system      (filename = os.path.join('/home/fernando/programs/pepdice/Examples/teste/RO1aiu/' , _file   )   )   

        system.ff_type = 'Calpha_model'
        system.import_Calpha_model_parameters(top       = None   ,   
                                              torsions  = os.path.join(PEPDICE_PARAMETER, 'amber/AMBER_rotamers.dat') )  

        print '%20s %20.10f %10s'%(_file, system.energy(log   = True,  AB_energy = True), rmsd[_file])
#--------------------------------------------------------------------------------------------------------------------------------------------
'''



'''

pdbs = [
        '/home/fernando/programs/pepdice/Examples/PDBs/alpha/1GAB/1GAB_A_AMBER.pdb']

system  = build_fragment_library_from_pdbs (
                                            molecule             = system ,
                                            frag_size            = 3      ,
                                            number_of_fragments  = 100    ,
                                            pdblist              = pdbs   ,
                                            )


import pickle

pickle.dump(system.fragments, open( "1gab_fragments.p", "wb" ) )















import pickle
system.fragments = pickle.load( open( "1gab_fragments.p", "rb" ) )

i =1
monte_carlo(molecule           = system      ,
            #random             = random      ,
            temperature        = 10000          ,
            Kb                 = 1           ,
            angle_range        = 5           ,
            nSteps             = 5000        ,
            fragment_rate      = 1.0         , #between 0  and 1
            fragment_sidechain = True        ,
            PhiPsi_rate        = 1.0         ,
            trajectory         = TRAJECTORY+str(i),
            pn                 = 1                       )

'''
