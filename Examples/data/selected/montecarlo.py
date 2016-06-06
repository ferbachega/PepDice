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

import random
random.seed(10)

pdbs = {
        '2n1p':'HSVSHARPRWFWFSLLLLAAGVGIYLLPNR'                                                 ,
        '2n2r':'QKLCQRPSGTWSGVCGNNNACKNQCIRLEKARHGSCNYVFPAHKCICYFPC'                            ,
        '2n2s':'SCGSECAPEPDCWGCCLVQCAPSICAGWCGGS'                                               ,
        '2n52':'VFSQGGQVDCGEFQDTKVYCTRESNPHCGSDGQTYGNKCAFCKAIVKSGGKISLKHPGKC'                   ,
        '2n5j':'GPHMTSELQMKVDFFRKLGYSSSEIHSVLQKLGVQADTNTVLGELVKHG'                              ,
        '2n5u':'MASWSHPQFEKIEGRMDVGQKVRVCRIRDRVAQDIIQKLGQVGQITGFKMTDGSGVGVIVTFDDRSSTWFFEDEVEVVG',
        '2n7f':'RDCQEKWEYCIVPILGFVYCCPGLICGPFVCV'                                               ,
        '2n7i':'GSFTMNDTTVWISVAVLSAVICLIIVWAVALKGYSMV'                                          ,
        '2n9c':'MTEYKLVVVGAGGVGKSHVW'                                                           ,
        '2nat':'KYEITTIHNLFRKLTHRLFRRNFGYTLR'                                                   ,
        '2nav':'HGEGTFTSDCSKQCEEGIGHKYPFCHCR'                                                   ,
        '2nbd':'MQIFVKTLTGKTITLEVEPSDTIENAKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG'   ,
        }

for pdb in pdbs:
    #-------------------------------------------------------------------#
    system = Molecule(name = pdb)                                       #
    system.build_peptide_from_sequence (                                #
                                        sequence    = pdbs[pdb]  ,      #
                                        _type       = 'amber'    ,      #
                                        force_field = 'ff03ua'   ,      #
                                        overwrite   = True       ,      #
                                        )                               #
    #-------------------------------------------------------------------#

    #-------------------------------------------------------------------#
    print system.energy()                                               #
    #'''                                                                #
    minimize(molecule = system,                                         #
             imin     = 1     ,                                         #
             maxcyc   = 100  ,                                         #
             ncyc     = 50   ,                                         #
             cut      = 99    ,                                         #
             rgbmax   = 10    ,                                         #
             igb      = 1     ,                                         #
             ntb      = 0     ,                                         #
             ntpr     = 1     ,                                         #
             ntr      = 0     )                                         #
    #'''                                                                #
    print system.energy()                                               #
    #-------------------------------------------------------------------#

    #-------------------------------------------------------------------------------
    #try:
    import pickle
    system.fragments = pickle.load( open( pdb + '/' +pdb+"_fragments_raw_5.p", "rb" ))
    ##-------------------------------------------------------------------------------
    monte_carlo(molecule           = system                                    ,
                random             = random                                     ,
                temperature        = 500                                       ,
                Kb                 = 1                                         ,
                angle_range        = 60                                        ,
                nSteps             = 200                                       ,
                fragment_rate      = 1.0                                       ,
                fragment_sidechain = False                                     ,
                log_frequence      = 10                                        ,
                PhiPsi_rate        = 0.0                                       ,
                trajectory         = pdb + '/' +pdb+'MC_raw_trajectory_5'      ,
                pn                 = 1                                         )
    #except:
    #    print 'fail ', pdb 

    system.fragments = pickle.load( open( pdb + '/' +pdb+"_fragments_min_5.p", "rb" ))
    monte_carlo(molecule           = system                                    ,
                random             = random                                     ,
                temperature        = 500                                        ,
                Kb                 = 1                                         ,
                angle_range        = 60                                        ,
                nSteps             = 200                                       ,
                fragment_rate      = 1.0                                       ,
                fragment_sidechain = False                                     ,
                log_frequence      = 10                                        ,
                PhiPsi_rate        = 0.0                                       ,
                trajectory         = pdb + '/' +pdb+'MC_min_trajectory_5'      ,
                pn                 = 1                                         )


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

'''
monte_carlo(molecule           = system                     ,
            temperature        = 1000                       ,
            Kb                 = 1                          ,  #0.0019872041               ,
            angle_range        = 60                         ,
            nSteps             = 10000                      ,
            fragment_rate      = 0.2                        ,
            fragment_sidechain = False                      ,
            log_frequence      = 10                         ,
            PhiPsi_rate        = 1.0                        ,
            trajectory         = 'MonteCarlo_trajectory.xyz',
            pn                 = 1                          )
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


