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
system.load_PDB_to_system      (filename = os.path.join(PEPDICE_EXAMPLES , 'data/1gab_amber03ua.pdb'   )   )   
system.import_AMBER_parameters(top       = os.path.join(PEPDICE_EXAMPLES , 'data/1gab_amber03ua.prmtop')   ,   
                                torsions = os.path.join(PEPDICE_PARAMETER, 'amber/AMBER_rotamers.dat') )   

TRAJECTORY  = os.path.join(PEPDICE_EXAMPLES , 'outputs/1gab_amber_example04_refold.xyz')
#---------------------------------------------------------------------------------------------------------


#------------------------------------------Extending---------------------------------------------------
for i in range (0,len(system.residues)):
    phi_final_angle = set_phi_psi_dihedral( molecule=system, resi=i, bond='PHI',  angle = 180)
    psi_final_angle = set_phi_psi_dihedral( molecule=system, resi=i, bond='PSI',  angle = 180)
    ome_final_angle = set_phi_psi_dihedral( molecule=system, resi=i, bond='OMEGA',angle = 180)
    save_XYZ_to_file(system,TRAJECTORY)

save_PDB_to_file(system,  os.path.join(PEPDICE_EXAMPLES , 'outputs/1gab_amber_example04_ff03ua_extended.pdb'))
#------------------------------------------------------------------------------------------------------

phi_and_psi_table = [
                     [ 9999.000,   92.403,  179.912],
                     [ -137.098,   46.270,  179.512],
                     [  -64.725,  -26.484,  179.271],
                     [  -50.710,  -23.944, -179.954],
                     [  -44.522,  -62.936, -178.637],
                     [  -30.471,  -38.510,  179.220],
                     [   69.944, -122.775, -179.814],
                     [   66.397,   36.128, -179.953],
                     [ -131.277,  -49.205, -179.895],
                     [  -80.433,  -15.370, -179.440],
                     [  -63.444,  -32.599,  179.988],
                     [  -92.361,  -44.116, -179.837],
                     [  -65.936,  -18.711, -179.912],
                     [  -99.000,  -43.829,  179.915],
                     [  -62.833,  -29.041, -178.916],
                     [  -68.038,  -15.552, -179.691],
                     [ -116.985,  -14.412, -179.951],
                     [  -90.210,  -32.996,  179.889],
                     [  -56.841,  -19.542,  179.491],
                     [  -86.953,  -31.221,  179.731],
                     [  -54.652,  -20.247,  179.999],
                     [  116.138,  -39.157, -179.967],
                     [  -54.026,  107.891, -179.885],
                     [ -115.203,    7.509,  179.937],
                     [  -59.946,  152.976, -179.882],
                     [ -100.458,   -9.132, -179.733],
                     [  -74.831,  -24.447,  179.866],
                     [  -94.954,  -42.391,  179.880],
                     [  -62.808,  -45.875, -179.844],
                     [  -65.778,  -39.549,  179.981],
                     [  -71.600,  -10.728, -179.869],
                     [  -65.937,  -24.360, -179.959],
                     [  -80.512,  -37.613, -179.824],
                     [  -57.881,  -20.441,  179.962],
                     [  -47.135,  155.935, -179.967],
                     [  -90.940,  -55.302, -179.891],
                     [ -125.018, -175.787, -179.966],
                     [  -63.210,  -30.321,  179.456],
                     [  -72.141,  -21.698,  179.390],
                     [  -86.478,  -31.095,  179.108],
                     [  -72.722,  -56.678,  178.853],
                     [  -62.864,  -27.934,  179.915],
                     [  -66.875,  -46.745,  179.748],
                     [  -59.199,  -35.154,  179.494],
                     [  -61.998,  -40.481,  179.691],
                     [  -65.896,  -20.692,  179.978],
                     [  -75.074,  -36.312,  179.420],
                     [  -58.470,  -27.412,  179.189],
                     [  -56.074,  -39.001, -179.401],
                     [ -109.629,  -62.240, -179.687],
                     [  -41.234,  -37.809, -179.547],
                     [ -115.616,  129.068, -179.397],
                     [   80.415, 9999.000, 9999.000],
                     ]
#------------------------------------------- refolding ------------------------------------------------
for i in range (0,len(system.residues)):
    phi_final_angle = set_phi_psi_dihedral( molecule=system, resi=i, bond='PHI',angle = phi_and_psi_table[i][0]) #angle = -57 )
    psi_final_angle = set_phi_psi_dihedral( molecule=system, resi=i, bond='PSI',angle = phi_and_psi_table[i][1]) #angle = -47 )
    ome_final_angle = set_phi_psi_dihedral( molecule=system, resi=i, bond='OMEGA',angle = phi_and_psi_table[i][2]) #angle = -47 )
    print system.residues[i].name ,  phi_final_angle, psi_final_angle, ome_final_angle
    save_XYZ_to_file(system,TRAJECTORY)

save_PDB_to_file(system,  os.path.join(PEPDICE_EXAMPLES , 'outputs/1gab_amber_example04_ff03ua_refold.pdb'))
#------------------------------------------------------------------------------------------------------















