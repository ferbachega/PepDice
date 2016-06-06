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


pdbs = {
        '2n1p':'HSVSHARPRWFWFSLLLLAAGVGIYLLPNR'                                                 ,
        '2n2r':'QKLCQRPSGTWSGVCGNNNACKNQCIRLEKARHGSCNYVFPAHKCICYFPC'                            ,
        '2n2s':'SCGSECAPEPDCWGCCLVQCAPSICAGWCGGS'                                               ,
        '2n52':'VFSQGGQVDCGEFQDTKVYCTRESNPHCGSDGQTYGNKCAFCKAIVKSGGKISLKHPGKC'                   ,
        '2n5j':'GPHMTSELQMKVDFFRKLGYSSSEIHSVLQKLGVQADTNTVLGELVKHG'                              ,
        '2n5u':'MASWSHPQFEKIEGRMDVGQKVRVCRIRDRVAQDIIQKLGQVGQITGFKMTDGSGVGVIVTFDDRSSTWFFEDEVEVVG',
        #'2n7f':'RDCQEKWEYCIVPILGFVYCCPGLICGPFVCV'                                               ,
        '2n7i':'GSFTMNDTTVWISVAVLSAVICLIIVWAVALKGYSMV'                                          ,
        '2n9c':'MTEYKLVVVGAGGVGKSHVW'                                                           ,
        '2nat':'KYEITTIHNLFRKLTHRLFRRNFGYTLR'                                                   ,
        '2nav':'HGEGTFTSDCSKQCEEGIGHKYPFCHCR'                                                   ,
        '2nbd':'MQIFVKTLTGKTITLEVEPSDTIENAKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG'   ,
        }

for pdb in pdbs:
    #--------------------------------------------------------------------------------------------------------#
    system = Molecule()                                                                                      #
    system.load_PDB_to_system      (filename = pdb + '/' + pdb + '_1_ff03ua.pdb' )                           #
    system.import_AMBER_parameters(top       = pdb + '/' + pdb + '_1_ff03ua.top' ,                           #
                                    torsions = os.path.join(PEPDICE_PARAMETER, 'amber/AMBER_rotamers.dat'))  #
    #--------------------------------------------------------------------------------------------------------#

    
    pdbfiles = os.listdir(pdb)
    pdbs_raw = []
    pdbs_min = []
    
    #-------------------------------------------------------------------------------------------------------#
    for pdbfile in pdbfiles:                                                                                #
        if 'minimized' in pdbfile:                                                                          #
            pdbs_min.append(pdb+'/'+pdbfile)                                                                #
        else:                                                                                               #
            if '.pdb' in pdbfile:
                pdbs_raw.append(pdb+'/'+pdbfile)                                                                #
    #-------------------------------------------------------------------------------------------------------#
    
    '''
    for _file in pdbs_min:
        print _file
        arq = open(_file, 'r')
        arq = arq.readline()
        if 'nan' in arq:
            os.remove(_file)
            print 'removing ', _file
    '''
    #print pdbs_raw
    #print pdbs_min 
    
    import pickle
    #system.energy(log = True)
    #try:
    system  = build_fragment_library_from_pdbs (
                                                molecule             = system   ,
                                                frag_size            = 5        ,
                                                number_of_fragments  = 100      ,
                                                pdblist              = pdbs_raw ,
                                               )
    pickle.dump(system.fragments, open( pdb + '/' +pdb+"_fragments_raw_5.p", "wb" ) )
    #except:
    #    print 'failed: ', pdb + '/' +pdb+"_fragments_raw.p"
    

    #try:
    system  = build_fragment_library_from_pdbs (
                                                molecule             = system   ,
                                                frag_size            = 5        ,
                                                number_of_fragments  = 100      ,
                                                pdblist              = pdbs_min ,
                                            )
    pickle.dump(system.fragments, open( pdb + '/' +pdb+"_fragments_min_5.p", "wb" ) )
    #except:
    #    print 'failed: ', pdb + '/' +pdb+"_fragments_min.p"

