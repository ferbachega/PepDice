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
system.load_PDB_to_system      (filename = os.path.join(PEPDICE_EXAMPLES , 'outputs/example04_1gab_refold.pdb'   )   )  
system.import_AMBER_parameters(top       = os.path.join(PEPDICE_EXAMPLES , 'data/alpha/1GAB/1gab_ff03ua_AMBER_extended.prmtop')   , 
                                torsions = os.path.join(PEPDICE_PARAMETER, 'amber/AMBER_rotamers.dat') )   

#TRAJECTORY  = os.path.join(PEPDICE_EXAMPLES , 'outputs/1gab_amber_example04_refold.xyz')                  /home/farminf/Programas/pepdice/Examples/outputs/example04_1gab_extended.pdb
#--------------------------------------------------------------------------------------------------------- /home/farminf/Programas/pepdice/Examples/outputs/example04_1gab_refold.pdb
# /home/farminf/Programas/pepdice/Examples/outputs/1gab_ff03ua_AMBER_extended_output.pdb
# /home/farminf/Programas/pepdice/Examples/data/alpha/1GAB/1gab_0002.pdb_UFF_min

#------------------------------------------Extending---------------------------------------------------
#for i in range (0,len(system.residues)):
#    phi_final_angle = set_phi_psi_dihedral( molecule=system, resi=i, bond='PHI',  angle = 180)
#    psi_final_angle = set_phi_psi_dihedral( molecule=system, resi=i, bond='PSI',  angle = 180)
#    ome_final_angle = set_phi_psi_dihedral( molecule=system, resi=i, bond='OMEGA',angle = 180)
#    save_XYZ_to_file(system,TRAJECTORY)
#
#save_PDB_to_file(system,  os.path.join(PEPDICE_EXAMPLES , 'outputs/1gab_amber_example04_ff03ua_extended.pdb'))
#------------------------------------------------------------------------------------------------------

torsions = []
log = True
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


phi_and_psi_table = torsions
system.load_PDB_to_system      (filename = os.path.join(PEPDICE_EXAMPLES , 'outputs/example04_1gab_extended.pdb'   )   )   
for i in range (0,len(system.residues)):
    phi_final_angle = set_phi_psi_dihedral( molecule=system, resi=i, bond='PHI',angle = phi_and_psi_table[i][0]) #angle = -57 )
    psi_final_angle = set_phi_psi_dihedral( molecule=system, resi=i, bond='PSI',angle = phi_and_psi_table[i][1]) #angle = -47 )
    ome_final_angle = set_phi_psi_dihedral( molecule=system, resi=i, bond='OMEGA',angle = phi_and_psi_table[i][2]) #angle = -47 )
    print system.residues[i].name ,  phi_final_angle, psi_final_angle, ome_final_angle
save_PDB_to_file(system,  os.path.join(PEPDICE_EXAMPLES , 'outputs/example04.1_1gab_refolded.pdb'))
#------------------------------------------------------------------------------------------------------















