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


#-------------------------------------------------------------------------------
system = Molecule() 
system.load_PDB_to_system      (filename = os.path.join(PEPDICE_EXAMPLES , 'data/alpha/1GAB/1gab_ff03ua_AMBER_folded.pdb'   )   )  
system.import_AMBER_parameters(top       = os.path.join(PEPDICE_EXAMPLES , 'data/alpha/1GAB/1gab_ff03ua_AMBER_extended.prmtop')   ,  
                                torsions = os.path.join(PEPDICE_PARAMETER, 'amber/AMBER_rotamers.dat') )   
#TRAJECTORY  = '/home/farminf/Programas/PepDice/Examples/1GAB/amber/1gab_amber_side_chain_rand.xyz' /home/farminf/Programas/pepdice/Examples/data/alpha/1GAB/1gab_ff03ua_AMBR_folded.pdb
#-------------------------------------------------------------------------------

#print system.energy( log = True)


def test_computeTorsions (system = None, log =False):
    """ Function doc """
    
    if log:
        print '''
----------------------------------------------------------
|                   TESTING TORSIONS                     |
----------------------------------------------------------
----------------------------------------------------------
| RESIDUE       PHI            PSI             OMEGA     |
----------------------------------------------------------
'''     
    torsions = []
    for i in range (0,len(system.residues)):
        phi_final_angle = computePhiPsi (molecule=system, resi=i, bond='PHI')
        psi_final_angle = computePhiPsi (molecule=system, resi=i, bond='PSI')
        ome_final_angle = computePhiPsi (molecule=system, resi=i, bond='OMEGA')

        if phi_final_angle == None:
            phi_final_angle = 0

        if psi_final_angle == None:
            psi_final_angle = 0
            
        if ome_final_angle == None:
            ome_final_angle = 0
        
        if log:
            print "%s  %15.5f  %15.5f  %15.5f" %(system.residues[i].name , phi_final_angle, psi_final_angle, ome_final_angle)
        
        
        torsions.append([phi_final_angle, 
                         psi_final_angle, 
                         ome_final_angle])
        

    if log:
        print '''
----------------------------------------------------------
\n\n''' 
    return torsions




torsions  = test_computeTorsions   (system = system, log = False)
pprint(torsions)
