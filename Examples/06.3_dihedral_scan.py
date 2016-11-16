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
from MonteCarlo  import monte_carlo#,  monte_carlo_dic, MC_replica_exchange, run_MC_replica_exchange #
from XYZFiles    import save_XYZ_to_file             #
from CRDFiles    import load_CRD_from_file           #
from AATorsions  import ROTAMER_LIST                 #
                                                     #
from RMSD import compute_RMSD                        #
                                                     #
from Test        import *                            #
                                                     #
from GeometryOptimization import minimize            #
from Energy import save_PDB_to_file                  #
                                                     #
#----------------------------------------------------#
from SideChainRefine import optimize_side_chain
#----------------------------------------------------#
import random
random.seed(1234)


#-------------------------------------------------------------------------------
PEPDICE = os.environ.get('PEPDICE')
PEPDICE_EXAMPLES = os.path.join(PEPDICE, 'Examples')
PEPDICE_PARAMETER= os.path.join(PEPDICE, 'Parameters')
#-------------------------------------------------------------------------------
#/home/farminf/Programas/PepDice/Examples/outputs/1gab_amber_example04_extended.pdb

system = Molecule()
name = 'poliALa_from_sequence'
TRAJECTORY = name+'_traj'
system.name = name
system.build_peptide_from_sequence (
                                     sequence    = 'AAA',
                                     _type       = 'amber'       ,
                                     force_field = 'ff03ua.labio',
                                     overwrite   = True          ,
                                     )
system.set_energy_model('RAW')

                                                    #AAAAAAAA
system.import_SS_restraints_from_string (      ss = 'HHH', 
                                            w_ss  = '090', log= True)
print system.compute_SS_energy(log = True)

system.energy_components['SS_RESTRAINT']  [1]= 10.0 
system.energy_components['SIZE'          ][1]= 0.0
system.energy_components['CONSTANT'      ][1]= 0.0
system.energy_components['ANGLE'         ][1]= 0.0
system.energy_components['BOND'          ][1]= 0.0
system.energy_components['DIHED'         ][1]= 0.0
system.energy_components['EEL'           ][1]= 0.0
system.energy_components['EELEC'         ][1]= 0.0
system.energy_components['EGB'           ][1]= 0.0
system.energy_components['EKtot'         ][1]= 0.0
system.energy_components['EPtot'         ][1]= 0.0
system.energy_components['ESURF'         ][1]= 0.0
system.energy_components['Etot'          ][1]= 0.0
system.energy_components['NB'            ][1]= 0.0
system.energy_components['VDWAALS'       ][1]= 0.0
system.energy_components['CONTACT'       ][1]= 0.0
system.energy_components['AB_ENERGY'     ][1]= 0.0
system.energy_components['R_GYRATION'    ][1]= 0.0




'''
for i in range (-180, 180 ):
    try:
        #rotate_backbone(molecule=system, resi=1, bond='PHI'  , theta= 0.1)
        
        set_phi_psi_dihedral (molecule = system,
                              resi  = 1     ,
                              bond  = 'PHI'    ,
                              angle = i)
        
        
        phi_final_angle = computePhiPsi (molecule=system, resi=1, bond='PHI')
        energy = system.energy()
        print 'phi: ',i,  phi_final_angle, energy
    except:
        print 'impossible to rotate PHI'
'''

set_phi_psi_dihedral (molecule = system,
                              resi  = 1     ,
                              bond  = 'PHI'    ,
                              angle = 179)


for i in range (-180, 180 ):
    try:
        #rotate_backbone(molecule=system, resi=1, bond='PHI'  , theta= 0.1)
        
        set_phi_psi_dihedral (molecule = system,
                              resi  = 1     ,
                              bond  = 'PSI'    ,
                              angle = i)
        
        
        psi_final_angle = computePhiPsi (molecule=system, resi=1, bond='PSI')
        energy = system.energy()
        print 'psi: ',i,  psi_final_angle, energy
    
    except:
        print 'impossible to rotate PSI'
#'''





system.Status()


i = system.energy()
'''
minimize(
        molecule=system,
        imin  =   1,
        maxcyc=1000,
        ncyc  = 100,
        cut   =  10,
        rgbmax= 999,
        igb   =   1,
        ntb   =   0,
        ntpr  = 100,
        ntr   =   0)
'''
save_PDB_to_file(system, '03_estendida_A_AMBER_minimized.pdb')

f = system.energy()

print 'delta= ', f-i 
#system.load_PDB_to_system      (filename = os.path.join(PEPDICE_EXAMPLES , '03_estendida_A_AMBER_minimized.pdb'   )   )   
#system.energy(log = True)


import pickle

system.fragments = pickle.load( open( os.path.join(PEPDICE_EXAMPLES,'1gab_fragments5.p'), "rb" ) )




fragments = system.fragments
#pprint (fragments)
#print len(fragments)
#print len(fragments[0])
#n = 0 

#for resi in fragments:
#    k = 0
#    if resi == []:
#        print n, resi
#    for frag in resi: 
#        print 'Position: ',n , 'fragment index: ',k, 'Number of fragments : ',len(resi), 'fragment size : ', len(frag)
#        k += 1
#    n += 1
#

system.energy(log = True)

#monte_carlo    (molecule            = system                 ,
#                random              = random                 ,
#                
#                initial_temperature = 100                    ,
#                final_temperature   = 0.1                    ,
#                gamma               = False                  ,
#                
#                
#                Kb                  = 1                      , # 0.0083144621               ,
#                angle_range         = 45                     ,
#                fragment_rate       = 0.0                    , #between 0  and 1
#                fragment_sidechain  = False                  ,
#                PhiPsi_rate         = 1.0                    ,
#                                    
#                
#                simulated_annealing = 'exp'                  , # exp, linear
#                cycle_size          = 100                    ,
#                number_of_cycles    = 10                     ,
#                
#                log_frequence       = 1                      ,
#                trajectory          = TRAJECTORY   ,
#                pn                  = 1                      ,
#                log                 = True                   ,)
#
