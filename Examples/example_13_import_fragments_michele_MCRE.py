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
from MonteCarlo  import monte_carlo,  monte_carlo_dic, MC_replica_exchange, run_MC_replica_exchange #
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

#-------------------------------------------------------------------------------
PEPDICE = os.environ.get('PEPDICE')
PEPDICE_EXAMPLES = os.path.join(PEPDICE, 'Examples')
PEPDICE_PARAMETER= os.path.join(PEPDICE, 'Parameters')
#-------------------------------------------------------------------------------
#/home/farminf/Programas/PepDice/Examples/outputs/1gab_amber_example04_extended.pdb


#-----------------------------------------------------------------------------------------------------------------------------------#
system = Molecule()                                                                                                                 #
#system.load_PDB_to_system      (filename = os.path.join(PEPDICE_EXAMPLES , 'data/alpha/1GAB/1gab_ff03ua_AMBER_extended.pdb'   )   ) #
system.load_PDB_to_system      (filename = os.path.join(PEPDICE_EXAMPLES , 'data/example03_geometry_opt_1gab_extended.pdb'   )   ) #
system.import_AMBER_parameters(top       = os.path.join(PEPDICE_EXAMPLES , 'data/alpha/1GAB/1gab_ff03ua_AMBER_extended.prmtop')   , #
                                torsions = os.path.join(PEPDICE_PARAMETER, 'amber/AMBER_rotamers.dat') )                            #

TRAJECTORY = os.path.join(PEPDICE_EXAMPLES , 'outputs/example10_n'   )   
try:                                                                                                                                #
    os.remove(TRAJECTORY)                                                                                                           #
except:                                                                                                                             #
    pass                                                                                                                            #
#-----------------------------------------------------------------------------------------------------------------------------------#
#system.energy(log = True)

    

import pickle
system.fragments = pickle.load( open( "/home/farminf/Programas/pepdice/Examples/data/alpha/1GAB/template_library_1gab_5.pkl", "rb" ) )
#system.fragments = pickle.load( open( "1gab_fragments.p", "rb" ) )

fragments = system.fragments
#pprint (fragments)
print len(fragments)
print len(fragments[0])
#n = 0 
#for resi in fragments:
#    k = 0
#    if resi == []:
#        print n, resi
#    for frag in resi: 
#        print 'Position: ',n , 'fragment index: ',k, 'Number of fragments : ',len(resi), 'fragment size : ', len(frag)
#        k += 1
#    n += 1


system.bond      = 1.0
system.angle     = 1.0
system.dihed     = 1.0
system.imprp     = 1.0
system.elect     = 1.0
system.vdw       = 1.0
system.boundary  = 1.0
system.esurf     = 1.0
system.egb       = 1.0


# usando o arquivo geran
#system.load_PDB_to_system      (filename = '/home/farminf/Documents/1GAB/1gab_pymol_refmac2.pdb') 
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
                        nSteps             = 500                 ,
                        nExchanges         = 5                   ,
                        
                        fragment_rate      = 0.3                 ,
                        log                = False               ,
                        #filelog            = 'MC_1GAB_replica_'
                        fragment_sidechain = True               ,
                        )
