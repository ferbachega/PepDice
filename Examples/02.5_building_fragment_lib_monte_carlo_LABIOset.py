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
from MonteCarlo  import monte_carlo, insert_fragment, insert_fragment_from_dic#,  monte_carlo_dic, MC_replica_exchange, run_MC_replica_exchange #
from XYZFiles    import save_XYZ_to_file                                                            #
from CRDFiles    import load_CRD_from_file                                                          #
from AATorsions  import ROTAMER_LIST                                                                #
                                                                                                    #
from RMSD import compute_RMSD                                                                       #
                                                                                                    #
from Test        import *                                                                           #
                                                                                                    #
from GeometryOptimization import minimize                                                           #
                                                                                                    #
from random import randint                                                                          #
                                                                                                    #
from Energy import save_PDB_to_file                                                                 #
                                                                                                    #
#---------------------------------------------------------------------------------------------------#
from SideChainRefine import optimize_side_chain
#----------------------------------------------------#
from Fragments import build_fragment_library_from_pdbs

import random as random 
random.seed(12345) 


#-------------------------------------------------------------------------------
PEPDICE = os.environ.get('PEPDICE')
PEPDICE_EXAMPLES = os.path.join(PEPDICE, 'Examples')
PEPDICE_PARAMETER= os.path.join(PEPDICE, 'Parameters')
#-------------------------------------------------------------------------------




# building a contact map - required for Contact model calculations
#-------------------------------------------------------------------------------
#from CMAP import CMAP
#cmap = CMAP(pdb = os.path.join(PEPDICE_EXAMPLES , 'LABIO_set/1I6C/1I6C_A_AMBER_minimized.pdb'), cutoff = 6.5, log = True)
#-------------------------------------------------------------------------------


def get_sequence_from_segFile (seqFile = None):
    """ Function doc """
    seqFile  =  open(seqFile, 'r')
    sequence =  False
    
    for line in seqFile:
        line2 = line.split()
        
        if len(line2) = 2:
            sequence = line2[1]



pdbs = [
        '1BDD',
        '1CSK',
        '1DV0',
        '1E0L',
        '1E0N',
        '1EI0',
        '1ERY',
        '1FME',
        '1GAB',
        '1I6C',
        '1L2Y',
        '1PRB',
        '1PSV',
        '1RES',
        '1RIJ',
        '1VII',
        '1WY3',
        '1YRF',
        '2HBA',
        '2HEP',
        '2JOF',
        '2WXC',
       ]







for pdb in pdbs:


    PDB        = pdb
    sequence   = get_sequence_from_segFile(seqFile = os.path.join(PEPDICE_EXAMPLES , '') )
    TRAJECTORY = PDB+'trajectory'
    # creating a new system 
    system = Molecule()
    system.name = PDB+'-from_sequence' 
    system.build_peptide_from_sequence (
                                         sequence    = sequence,
                                         _type       = 'amber'       ,
                                         force_field = 'ff03ua.labio',
                                         overwrite   = True          ,
                                         )
    system.set_energy_model('RAW')
    system.load_PDB_to_system       (filename =  os.path.join(PEPDICE_EXAMPLES , 'new_method_example02.1_fragments_extended_min.pdb'))



    '''
    minimize(molecule = system,
                   imin  = 1          ,
                   maxcyc= 1000       ,
                   ncyc  = 100        ,
                   cut   = 10         ,
                   rgbmax= 999        ,
                   igb   = 1          ,
                   ntb   = 0          ,
                   ntpr  = 100        ,
                   ntr   = 0          )
                   #restraintmask = ':1-50 & @CA,N,C,O=', 
                   #restraint_wt  =  50.0               )
    save_PDB_to_file(system,  os.path.join(PEPDICE_EXAMPLES , 'new_method_example02.1_fragments_extended_min.pdb'))


    system.load_PDB_to_system       (filename = os.path.join( PEPDICE_EXAMPLES , 'LABIO_set/'+PDB+'/'+PDB+'_A_AMBER_minimized.pdb'     ))
    '''


    decoys = 


    pdbs          = [
                        os.path.join( PEPDICE_EXAMPLES , 'LABIO_set/'+PDB+'/'+PDB+'_A_AMBER_minimized.pdb'     ),
                        os.path.join( PEPDICE_EXAMPLES , 'LABIO_set/'+PDB+'/decoy0_1_A_AMBER_minimized.pdb' ),
                        os.path.join( PEPDICE_EXAMPLES , 'LABIO_set/'+PDB+'/decoy1_1_A_AMBER_minimized.pdb' ),
                        #os.path.join( PEPDICE_EXAMPLES , 'LABIO_set/'+PDB+'/decoy1_2_A_AMBER_minimized.pdb' ),
                        #os.path.join( PEPDICE_EXAMPLES , 'LABIO_set/'+PDB+'/decoy1_3_A_AMBER_minimized.pdb' ),
                        #os.path.join( PEPDICE_EXAMPLES , 'LABIO_set/'+PDB+'/decoy1_4_A_AMBER_minimized.pdb' ),
                        #os.path.join( PEPDICE_EXAMPLES , 'LABIO_set/'+PDB+'/decoy1_5_A_AMBER_minimized.pdb' ),
                        #os.path.join( PEPDICE_EXAMPLES , 'LABIO_set/'+PDB+'/decoy1_6_A_AMBER_minimized.pdb' ),
                        #os.path.join( PEPDICE_EXAMPLES , 'LABIO_set/'+PDB+'/decoy1_7_A_AMBER_minimized.pdb' ),
                        #os.path.join( PEPDICE_EXAMPLES , 'LABIO_set/'+PDB+'/decoy1_8_A_AMBER_minimized.pdb' ),
                        os.path.join( PEPDICE_EXAMPLES , 'LABIO_set/'+PDB+'/decoy1_9_A_AMBER_minimized.pdb' ),
                        os.path.join( PEPDICE_EXAMPLES , 'LABIO_set/'+PDB+'/decoy1_10_A_AMBER_minimized.pdb'),
                        os.path.join( PEPDICE_EXAMPLES , 'LABIO_set/'+PDB+'/decoy1_11_A_AMBER_minimized.pdb'),
                        #os.path.join( PEPDICE_EXAMPLES , 'LABIO_set/'+PDB+'/decoy1_12_A_AMBER_minimized.pdb'),
                        #os.path.join( PEPDICE_EXAMPLES , 'LABIO_set/'+PDB+'/decoy1_13_A_AMBER_minimized.pdb'),
                        #os.path.join( PEPDICE_EXAMPLES , 'LABIO_set/'+PDB+'/decoy1_14_A_AMBER_minimized.pdb'),
                        #os.path.join( PEPDICE_EXAMPLES , 'LABIO_set/'+PDB+'/decoy1_15_A_AMBER_minimized.pdb'),
                        #os.path.join( PEPDICE_EXAMPLES , 'LABIO_set/'+PDB+'/decoy1_16_A_AMBER_minimized.pdb'),
                        #os.path.join( PEPDICE_EXAMPLES , 'LABIO_set/'+PDB+'/decoy1_17_A_AMBER_minimized.pdb'),
                        #os.path.join( PEPDICE_EXAMPLES , 'LABIO_set/'+PDB+'/decoy1_23_A_AMBER_minimized.pdb'),
                        #os.path.join( PEPDICE_EXAMPLES , 'LABIO_set/'+PDB+'/decoy1_24_A_AMBER_minimized.pdb'),
                        #os.path.join( PEPDICE_EXAMPLES , 'LABIO_set/'+PDB+'/decoy1_30_A_AMBER_minimized.pdb'),
                        os.path.join( PEPDICE_EXAMPLES , 'LABIO_set/'+PDB+'/decoy1_31_A_AMBER_minimized.pdb'),
                        os.path.join( PEPDICE_EXAMPLES , 'LABIO_set/'+PDB+'/decoy2_14_A_AMBER_minimized.pdb'),
                        os.path.join( PEPDICE_EXAMPLES , 'LABIO_set/'+PDB+'/decoy2_35_A_AMBER_minimized.pdb')
                        ]



    from Fragments import *

    system  = build_fragment_library_from_pdbs (
                                                molecule             = system ,
                                                frag_size            = 7      ,
                                                number_of_fragments  = 20     ,
                                                pdblist              = pdbs   ,
                                                )


    system.Status()
    system.load_PDB_to_system       (filename =  os.path.join(PEPDICE_EXAMPLES , 'new_method_example02.1_fragments_extended_min.pdb'))


    '''
    monte_carlo    (molecule            = system                 ,
                    random              = random                 ,
                    
                    initial_temperature = 100                    ,
                    final_temperature   = 0.1                    ,
                    gamma               = False                  ,
                    
                    
                    Kb                  = 1                      , # 0.0083144621               ,
                    angle_range         = 1.0                   ,
                    fragment_rate       = 1.0                    , #between 0  and 1
                    fragment_sidechain  = True                   ,
                    PhiPsi_rate         = 1.0                    ,
                                        
                    
                    simulated_annealing = 'exp'                  , # exp, linear
                    cycle_size          = 100                    ,
                    number_of_cycles    = 10                     ,
                    
                    log_frequence       = 1                      ,
                    trajectory          = TRAJECTORY   ,
                    pn                  = 1                      ,
                    log                 = True                   ,)

    '''
