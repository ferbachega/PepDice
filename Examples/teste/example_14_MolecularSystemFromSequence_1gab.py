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
from Energy import save_PDB_to_file                                                                 #
                                                                                                    #
#---------------------------------------------------------------------------------------------------#
from SideChainRefine import optimize_side_chain                                                     #
#---------------------------------------------------------------------------------------------------#

#---------------------------------------------------------------------------------------------------
import random
random.seed(10)                                                                                     
#---------------------------------------------------------------------------------------------------                                                                                      #
                                                                                                    



#---------------------------------------------------------------------------------------------------
system = Molecule(name = '1GAB_poly_ala') 
system.build_peptide_from_sequence (
                                    sequence    = 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA',
                                    _type       = 'amber'    ,
                                    force_field = 'ff03ua'   ,
                                    overwrite   = True       ,
                                    )
#---------------------------------------------------------------------------------------------------#

'''
print 'name    : ',system.name    
print 'top     : ',system.top     
print 'psf     : ',system.psf     
print 'param   : ',system.param   
print 'ff_type : ',system.ff_type 
print 'torsions: ',system.torsions 
#'''


system.bond      = 1.0
system.angle     = 1.0
system.dihed     = 1.0
system.imprp     = 1.0
system.elect     = 1.0
system.vdw       = 1.0
system.boundary  = 1.0
system.esurf     = 1.0
system.egb       = 1.0



#-------------------------------------------------------------------------------
print system.energy()
#'''
minimize(molecule = system,
         imin     = 1     ,
         maxcyc   = 100   ,
         ncyc     = 100   ,
         cut      = 99    ,
         rgbmax   = 10    ,
         igb      = 1     ,
         ntb      = 0     ,
         ntpr     = 1     ,
         ntr      = 0     )
#'''
print system.energy()
#-------------------------------------------------------------------------------

import pickle
#system.fragments = pickle.load( open( "1gab_fragments_bachega.p", "rb" ) )
system.fragments = pickle.load( open( "1gab_fragments.p", "rb"))#/home/farminf/Programas/pepdice/Examples/data/alpha/1GAB/template_library_1gab_5.pkl", "rb" ) )


'''
run_MC_replica_exchange (
                        molecule           = system              ,
                        N_replicas         = 8                   , # >= number of CPUs
                        CPUs               = 8                   ,
                        min_temp           = 100                 ,
                        max_temp           = 1000                ,
                        PhiPsi_rate        = 1.0                 ,
                        max_angle_range    = 5                   ,
                        trajectory         = 'MC_1GAB_replica_'  ,
                        Kb                 = 0.0019872041        ,
                        log_frequence      = 10                  ,
                        nSteps             = 10000               ,
                        nExchanges         = 5                   ,

                        fragment_rate      = 0                 ,
                        log                = False               ,
                        #filelog            = 'MC_1GAB_replica_'
                        fragment_sidechain = True               ,
                        )
'''


#'''
monte_carlo(molecule           = system                        ,
            temperature        = 10                            ,
            Kb                 = 1                             ,  #0.0019872041               ,
            angle_range        = 60                            ,
            nSteps             = 10000                         ,
            fragment_rate      = 1.0                           ,
            fragment_sidechain = False                         ,
            log_frequence      = 10                            ,
            PhiPsi_rate        = 0.0                           ,
            trajectory         = system.name + '_MC_trajectory',
            pn                 = 1                             ,
            random             = random                        )
#'''


#TRAJECTORY  = os.path.join(PEPDICE_EXAMPLES , 'outputs/1gab_amber_example04_refold.xyz')
#---------------------------------------------------------------------------------------------------------


'''
def test_rotetaPhiPsiOMEGA (system = None, theta =  0.017453292519943295,  computeTorsion = False):
    """ Function doc """
    backup  = system 
    
    for j in range(0,10):
        system = backup
        for i in range (0,len(system.residues)):

            try:
                rotate_backbone(molecule=system, resi=i, bond='PHI'  , theta= theta*(randint(j*-10,j)))
            except:
                print 'impossible to rotate PHI'
            
            try:
                rotate_backbone(molecule=system, resi=i, bond='PSI'  , theta=theta*(randint(j*-10,j)))
            except:
                print 'impossible to rotate PSI'

            #try:
            #    rotate_backbone(molecule=system, resi=i, bond='OMEGA', theta=theta*(randint(j*-1,j)))
            #except:
            #    print 'impossible to rotate OMEGA'
                
                
        save_PDB_to_file(system, 'example05_polyAla_rotetaBackbone_'+str(j)+'_rand.pdb')
        #test_computeTorsions (system = system)

test_rotetaPhiPsiOMEGA (system = system)
#'''


