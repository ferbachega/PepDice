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


#---------------------------------------------------------------------------------------------------------
system = Molecule(name = '1gab') 
system.build_peptide_from_sequence (
                                     sequence    = 'AAAAAAAAA',
                                     _type       = 'amber'    ,
                                     force_field = 'ff03ua'   ,
                                     overwrite   = True       ,
                                     )




system.bond      = 1.0
system.angle     = 1.0
system.dihed     = 1.0
system.imprp     = 1.0
system.elect     = 1.0
system.vdw       = 1.0
system.boundary  = 1.0
system.esurf     = 1.0
system.egb       = 1.0



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

                        fragment_rate      = 1.0                 ,
                        log                = False               ,
                        #filelog            = 'MC_1GAB_replica_'
                        fragment_sidechain = True               ,
                        )








