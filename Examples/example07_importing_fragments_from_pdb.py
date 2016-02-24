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

TRAJECTORY  = os.path.join(PEPDICE_EXAMPLES , 'outputs/1gab_amber_example04_refold.xyz')
#---------------------------------------------------------------------------------------------------------

def import_fragments_from_pdb (molecule = None, residues = []):
    """ Function doc """
    fragment = {}
    
    for i in residues:
        fragment[i] = {}
        fragment[i]['PHI']   = computePhiPsi (molecule=system, resi=i, bond='PHI')
        fragment[i]['PSI']   = computePhiPsi (molecule=system, resi=i, bond='PSI')
        fragment[i]['OMEGA'] = computePhiPsi (molecule=system, resi=i, bond='OMEGA')
      
        for chi in ["CHI1","CHI2","CHI3","CHI4","CHI5"]:
            fragment[i][chi] = computeCHI (molecule= molecule, resi=i, bond=chi)
            #print i, system.residues[i].name, chi,computeCHI (molecule=system, resi=i, bond=chi)
    
    
    return fragment 


fragments = []
for i in range(0, 10):
    resi = random.randint(0,40)
    fragment = import_fragments_from_pdb(molecule = system, residues= range(resi, resi+7))
    #pprint (fragment)
    fragments.append(fragment)
pprint(fragments)

'''
#------------------------------------------Extending---------------------------------------------------
for i in range (0,len(system.residues)):
    for chi in ["CHI1","CHI2","CHI3","CHI4","CHI5"]:
        print i, system.residues[i].name, chi,computeCHI (molecule=system, resi=i, bond=chi)
    #phi_final_angle = set_phi_psi_dihedral( molecule=system, resi=i, bond='PHI',  angle = 180)
    #psi_final_angle = set_phi_psi_dihedral( molecule=system, resi=i, bond='PSI',  angle = 180)
    #ome_final_angle = set_phi_psi_dihedral( molecule=system, resi=i, bond='OMEGA',angle = 180)
    #save_XYZ_to_file(system,TRAJECTORY)

#save_PDB_to_file(system,  os.path.join(PEPDICE_EXAMPLES , 'outputs/1gab_amber_example04_extended.pdb'))
#------------------------------------------------------------------------------------------------------
'''













