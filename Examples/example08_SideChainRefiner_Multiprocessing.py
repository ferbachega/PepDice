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

#-------------------------------------------------------------------------------
PEPDICE = os.environ.get('PEPDICE')
PEPDICE_EXAMPLES = os.path.join(PEPDICE, 'Examples')
PEPDICE_PARAMETER= os.path.join(PEPDICE, 'Parameters')
#-------------------------------------------------------------------------------


#---------------------------------------------------------------------------------------------------------
system = Molecule() 
system.load_PDB_to_system      (filename = os.path.join(PEPDICE_EXAMPLES , 'data/1gab_amber.pdb'   )   )   
system.import_AMBER_parameters(top       = os.path.join(PEPDICE_EXAMPLES , 'data/1gab_amber.prmtop')   ,   
                                torsions = os.path.join(PEPDICE_PARAMETER, 'amber/AMBER_rotamers.dat') )   
TRAJECTORY  = os.path.join(PEPDICE_EXAMPLES , 'outputs/1gab_amber_example08_sidechain_refine.xyz')
#---------------------------------------------------------------------------------------------------------


system.energy(log = True)


#energy     = True
#trajectory = TRAJECTORY
#initial_coordinates = system.get_coordinates_from_system()



def optimize_pairs (system = None, i = None, j = None):
    """ Function doc """
    name         = system.residues[i].name
    res          = system.residues[i]
    
    if name in ['HSD', 'HSE', 'HDP', 'HIE', 'HID']:
        name = 'HIS'
    
    res_rotamers =  ROTAMER_LIST[name]
    energy_list = {}

    for key in res_rotamers:
        set_side_chain_rotamer(molecule=system, resi=i, rotamer=res_rotamers[key])
        results = optimize_side_chain(system = system, resi = j)
        keys    = results.keys()
        try:
            lower_E = min(keys)
            print i, j, lower_E
            save_XYZ_to_file(system, TRAJECTORY)
            system.import_coordinates_to_system(results[lower_E])
            energy_list[lower_E] =  system.get_coordinates_from_system()
        except:
            print i, j, 'fail'

    try:
        lista = energy_list.keys()
        lower_E  = min(lista)
        coord    = energy_list[lower_E]
        system.import_coordinates_to_system(coord)
        save_XYZ_to_file(system, TRAJECTORY)
    except:
        pass
        
    
for resi in range(0,len(system.residues)):
    optimize_pairs (system = system, i = resi, j = resi+1)
