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




#-----------------------------------------------------------------------------------------------------------------------------------#
system = Molecule()                                                                                                                 #
system.load_PDB_to_system      (filename = os.path.join(PEPDICE_EXAMPLES , 'data/alpha/1GAB/1gab_ff03ua_AMBER_folded.pdb'     )   ) #
system.import_AMBER_parameters(top       = os.path.join(PEPDICE_EXAMPLES , 'data/alpha/1GAB/1gab_ff03ua_AMBER_extended.prmtop')   , #
                                torsions = os.path.join(PEPDICE_PARAMETER, 'amber/AMBER_rotamers.dat') )                            #
try:                                                                                                                                #
    os.remove(TRAJECTORY)                                                                                                           #
except:                                                                                                                             #
    pass                                                                                                                            #
#-----------------------------------------------------------------------------------------------------------------------------------#
system.energy(log = True)


pdbs          = [
                 os.path.join(PEPDICE_EXAMPLES ,'data/alpha/1GAB/1gab_0001.pdb'),
                 os.path.join(PEPDICE_EXAMPLES ,'data/alpha/1GAB/1gab_0002.pdb'),
                 os.path.join(PEPDICE_EXAMPLES ,'data/alpha/1GAB/1gab_0003.pdb'),
                 os.path.join(PEPDICE_EXAMPLES ,'data/alpha/1GAB/1gab_0004.pdb'),
                 os.path.join(PEPDICE_EXAMPLES ,'data/alpha/1GAB/1gab_0005.pdb'),
                 os.path.join(PEPDICE_EXAMPLES ,'data/alpha/1GAB/1gab_0006.pdb'),
                 os.path.join(PEPDICE_EXAMPLES ,'data/alpha/1GAB/1gab_0007.pdb'),
                 os.path.join(PEPDICE_EXAMPLES ,'data/alpha/1GAB/1gab_0008.pdb'),
                 os.path.join(PEPDICE_EXAMPLES ,'data/alpha/1GAB/1gab_0009.pdb'),
                 os.path.join(PEPDICE_EXAMPLES ,'data/alpha/1GAB/1gab_0010.pdb'),
                 os.path.join(PEPDICE_EXAMPLES ,'data/alpha/1GAB/1gab_0011.pdb'),
                 os.path.join(PEPDICE_EXAMPLES ,'data/alpha/1GAB/1gab_0012.pdb'),
                 os.path.join(PEPDICE_EXAMPLES ,'data/alpha/1GAB/1gab_0013.pdb'),
                 os.path.join(PEPDICE_EXAMPLES ,'data/alpha/1GAB/1gab_0014.pdb'),
                 os.path.join(PEPDICE_EXAMPLES ,'data/alpha/1GAB/1gab_0015.pdb'),
                 os.path.join(PEPDICE_EXAMPLES ,'data/alpha/1GAB/1gab_0016.pdb'),
                 os.path.join(PEPDICE_EXAMPLES ,'data/alpha/1GAB/1gab_0017.pdb'),
                 os.path.join(PEPDICE_EXAMPLES ,'data/alpha/1GAB/1gab_0018.pdb'),
                 os.path.join(PEPDICE_EXAMPLES ,'data/alpha/1GAB/1gab_0019.pdb'),
                 os.path.join(PEPDICE_EXAMPLES ,'data/alpha/1GAB/1gab_0020.pdb')]

system  = build_fragment_library_from_pdbs (
                                            molecule             = system ,
                                            frag_size            = 3      ,
                                            number_of_fragments  = 100    ,
                                            pdblist              = pdbs   ,
                                            )


import pickle

pickle.dump(system.fragments, open( "1gab_fragments.p", "wb" ) )


