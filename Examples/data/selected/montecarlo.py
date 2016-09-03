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
        '1gab': ['TIDQWLLKNAKEDAIAELKKAGITSDFYFNAINKAKTVEEVNALKNEILKAHA'                          , '1gab/1gab_fragments_min_5.p'],
        #'2n1p':['HSVSHARPRWFWFSLLLLAAGVGIYLLPNR'                                                 , None],
        #'2n2r':['QKLCQRPSGTWSGVCGNNNACKNQCIRLEKARHGSCNYVFPAHKCICYFPC'                            , None],
        #'2n2s':['SCGSECAPEPDCWGCCLVQCAPSICAGWCGGS'                                               , None],
        #'2n52':['VFSQGGQVDCGEFQDTKVYCTRESNPHCGSDGQTYGNKCAFCKAIVKSGGKISLKHPGKC'                   , None],
        #'2n5j':['GPHMTSELQMKVDFFRKLGYSSSEIHSVLQKLGVQADTNTVLGELVKHG'                              , None],
        #'2n5u':['MASWSHPQFEKIEGRMDVGQKVRVCRIRDRVAQDIIQKLGQVGQITGFKMTDGSGVGVIVTFDDRSSTWFFEDEVEVVG', None],
        #'2n7f':['RDCQEKWEYCIVPILGFVYCCPGLICGPFVCV'                                               , None],
        #'2n7i':['GSFTMNDTTVWISVAVLSAVICLIIVWAVALKGYSMV'                                          , None],
        #'2n9c':['MTEYKLVVVGAGGVGKSHVW'                                                           , None],
        #'2nat':['KYEITTIHNLFRKLTHRLFRRNFGYTLR'                                                   , None],
        #'2nav':['HGEGTFTSDCSKQCEEGIGHKYPFCHCR'                                                   , None],
        #'2nbd':['MQIFVKTLTGKTITLEVEPSDTIENAKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG'   , None],
        }

for pdb in pdbs:
    #-------------------------------------------------------------------#
    system = Molecule(name = pdb)                                       #
    system.build_peptide_from_sequence (                                #
                                        sequence    = pdbs[pdb][0],     #
                                        _type       = 'amber'     ,     #
                                        force_field = 'ff03ua'    ,     #
                                        overwrite   = True        ,     #
                                        )                               #
    #-------------------------------------------------------------------#

    #-------------------------------------------------------------------#
    print system.energy()                                               #
    #'''                                                                #
    minimize(molecule = system,                                         #
             imin     = 1     ,                                         #
             maxcyc   = 1000  ,                                         #
             ncyc     = 500   ,                                         #
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
    import pickle
    #-------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------
    system.fragments = pickle.load( open( pdbs[pdb][1], "rb" ))
    monte_carlo(molecule           = system                                    ,
                random             = random                                    ,
                temperature        = 1000                                       ,
                Kb                 = 1                                         ,
                angle_range        = 60                                        ,
                nSteps             = 10000                                     ,
                fragment_rate      = 1.0                                       ,
                fragment_sidechain = False                                     ,
                log_frequence      = 10                                        ,
                PhiPsi_rate        = 0.0                                       ,
                trajectory         = pdb + '/' +pdb+'MC_min_trajectory_5'      ,
                pn                 = 1                                         )
    #-------------------------------------------------------------------------------
    
    
    
    ##-------------------------------------------------------------------------------
    #system.fragments = pickle.load( open( pdb + '/' +pdb+"_fragments_min_5.p", "rb" ))
    #monte_carlo(molecule           = system                                    ,
    #            random             = random                                     ,
    #            temperature        = 1000                                       ,
    #            Kb                 = 1                                         ,
    #            angle_range        = 60                                        ,
    #            nSteps             = 10000                                     ,
    #            fragment_rate      = 1.0                                       ,
    #            fragment_sidechain = False                                     ,
    #            log_frequence      = 10                                        ,
    #            PhiPsi_rate        = 0.0                                       ,
    #            trajectory         = pdb + '/' +pdb+'MC_min_trajectory_5'      ,
    #            pn                 = 1                                         )
    ##-------------------------------------------------------------------------------


