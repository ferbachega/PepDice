#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  SideChainRefiner_Multiprocessing.py
#  
#  Copyright 2016 farminf <farminf@farminf-3>
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

#---------------------------------------------------------------------------------------------------#
import os                                                                                           #
from pprint      import pprint                                                                      #
from Molecule    import Molecule                                                                    #
from Geometry    import *                                                                           #
from MonteCarlo  import monte_carlo#,  monte_carlo_dic, MC_replica_exchange, run_MC_replica_exchange #
from XYZFiles    import save_XYZ_to_file                                                            #
from CRDFiles    import load_CRD_from_file                                                          #
from AATorsions  import ROTAMER_LIST                                                                #
                                                                                                    #
from RMSD import compute_RMSD                                                                       #
                                                                                                    #
from Test        import *                                                                           #
                                                                                                    #
from GeometryOptimization import minimize                                                           #
                                                                                                    #
from random import randint                                                                          #
                                                                                                    #
from Energy import save_PDB_to_file                                                                 #
                                                                                                    #
#---------------------------------------------------------------------------------------------------#
from SideChainRefine import optimize_side_chain
#----------------------------------------------------#
from Fragments import build_fragment_library_from_pdbs



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


pdbs          = [
                    os.path.join( PEPDICE_EXAMPLES , 'LABIO_set/1I6C/1I6C_A_AMBER_minimized.pdb'     ),
                    os.path.join( PEPDICE_EXAMPLES , 'LABIO_set/1I6C/decoy0_1_A_AMBER_minimized.pdb' ),
                    os.path.join( PEPDICE_EXAMPLES , 'LABIO_set/1I6C/decoy1_1_A_AMBER_minimized.pdb' ),
                    os.path.join( PEPDICE_EXAMPLES , 'LABIO_set/1I6C/decoy1_2_A_AMBER_minimized.pdb' ),
                    os.path.join( PEPDICE_EXAMPLES , 'LABIO_set/1I6C/decoy1_3_A_AMBER_minimized.pdb' ),
                    os.path.join( PEPDICE_EXAMPLES , 'LABIO_set/1I6C/decoy1_4_A_AMBER_minimized.pdb' ),
                    os.path.join( PEPDICE_EXAMPLES , 'LABIO_set/1I6C/decoy1_5_A_AMBER_minimized.pdb' ),
                    os.path.join( PEPDICE_EXAMPLES , 'LABIO_set/1I6C/decoy1_6_A_AMBER_minimized.pdb' ),
                    os.path.join( PEPDICE_EXAMPLES , 'LABIO_set/1I6C/decoy1_7_A_AMBER_minimized.pdb' ),
                    os.path.join( PEPDICE_EXAMPLES , 'LABIO_set/1I6C/decoy1_8_A_AMBER_minimized.pdb' ),
                    os.path.join( PEPDICE_EXAMPLES , 'LABIO_set/1I6C/decoy1_9_A_AMBER_minimized.pdb' ),
                    os.path.join( PEPDICE_EXAMPLES , 'LABIO_set/1I6C/decoy1_10_A_AMBER_minimized.pdb'),
                    os.path.join( PEPDICE_EXAMPLES , 'LABIO_set/1I6C/decoy1_11_A_AMBER_minimized.pdb'),
                    os.path.join( PEPDICE_EXAMPLES , 'LABIO_set/1I6C/decoy1_12_A_AMBER_minimized.pdb'),
                    os.path.join( PEPDICE_EXAMPLES , 'LABIO_set/1I6C/decoy1_13_A_AMBER_minimized.pdb'),
                    os.path.join( PEPDICE_EXAMPLES , 'LABIO_set/1I6C/decoy1_14_A_AMBER_minimized.pdb'),
                    os.path.join( PEPDICE_EXAMPLES , 'LABIO_set/1I6C/decoy1_15_A_AMBER_minimized.pdb'),
                    os.path.join( PEPDICE_EXAMPLES , 'LABIO_set/1I6C/decoy1_16_A_AMBER_minimized.pdb'),
                    os.path.join( PEPDICE_EXAMPLES , 'LABIO_set/1I6C/decoy1_17_A_AMBER_minimized.pdb'),
                    os.path.join( PEPDICE_EXAMPLES , 'LABIO_set/1I6C/decoy1_23_A_AMBER_minimized.pdb'),
                    os.path.join( PEPDICE_EXAMPLES , 'LABIO_set/1I6C/decoy1_24_A_AMBER_minimized.pdb'),
                    os.path.join( PEPDICE_EXAMPLES , 'LABIO_set/1I6C/decoy1_30_A_AMBER_minimized.pdb'),
                    os.path.join( PEPDICE_EXAMPLES , 'LABIO_set/1I6C/decoy1_31_A_AMBER_minimized.pdb'),
                    os.path.join( PEPDICE_EXAMPLES , 'LABIO_set/1I6C/decoy2_14_A_AMBER_minimized.pdb'),
                    os.path.join( PEPDICE_EXAMPLES , 'LABIO_set/1I6C/decoy2_35_A_AMBER_minimized.pdb')]















n = 3
for i in range(0,3):
    
    system  = build_fragment_library_from_pdbs (
                                                molecule             = system ,
                                                frag_size            = n      ,
                                                number_of_fragments  = 30     ,
                                                pdblist              = pdbs   ,
                                                )


    import pickle

    pickle.dump(system.fragments, open( "1gab_fragments"+str(n)+".p", "wb" ) )
    n += 2

