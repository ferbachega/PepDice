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
from MonteCarlo  import monte_carlo, insert_fragment #,  monte_carlo_dic, MC_replica_exchange, run_MC_replica_exchange #
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

from Fragments import import_fragments_from_pdb


#-------------------------------------------------------------------------------
PEPDICE = os.environ.get('PEPDICE')
PEPDICE_EXAMPLES = os.path.join(PEPDICE, 'Examples')
PEPDICE_PARAMETER= os.path.join(PEPDICE, 'Parameters')
#-------------------------------------------------------------------------------







#-----------------------------------------------------------------------------------------------------------------------------------#
system = Molecule()
#system.load_PDB_to_system      (filename = os.path.join(PEPDICE_EXAMPLES , 'data/alpha/1GAB/1gab_ff03ua_AMBER_folded.pdb'     )   )                                                                                                              #
system.load_PDB_to_system      (filename = os.path.join(PEPDICE_EXAMPLES , 'data/alpha/1GAB/1gab_ff03ua_AMBER_extended.pdb'   )   ) #
system.import_AMBER_parameters(top       = os.path.join(PEPDICE_EXAMPLES , 'data/alpha/1GAB/1gab_ff03ua_AMBER_extended.prmtop')   , #
                                torsions = os.path.join(PEPDICE_PARAMETER, 'amber/AMBER_rotamers.dat') )                            #

TRAJECTORY = os.path.join(PEPDICE_EXAMPLES , 'outputs/example10_novo_n'   )
try:                                                                                                                                #
    os.remove(TRAJECTORY)                                                                                                           #
except OSError:                                                                                                                             #
    # File doesnt exist, so nothing to do here
    pass                                                                                                                      #
#-----------------------------------------------------------------------------------------------------------------------------------#

system.energy(log = True)



#import pickle
#system.fragments = pickle.load( open( "1gab_fragments.p", "rb" ) )
#---------------------------------------------------------------------------------
#pdb     = os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy0_1.pdb')
pdb     = os.path.join(PEPDICE_EXAMPLES , 'data/alpha/1GAB/1gab_ff03ua_AMBER_folded.pdb'     )  
system2 = Molecule()                                                                                                                 #
system2.load_PDB_to_system(filename = pdb)





fragment = import_fragments_from_pdb (molecule  = system2, 
                                      residues  = range(0,53), 
                                      mainchain = True, 
                                      sidechain = False)
#---------------------------------------------------------------------------------




save_PDB_to_file(system,os.path.join(PEPDICE_EXAMPLES, 'data/01_inicial.pdb'))


minimize(molecule = system,
               imin  = 1          ,
               maxcyc= 1000        ,
               ncyc  = 100        ,
               cut   = 10         ,
               rgbmax= 999        ,
               igb   = 1          ,
               ntb   = 0          ,
               ntpr  = 100        ,
               ntr   = 0          )
               #restraintmask = ':1-50 & @CA,N,C,O=', 
               #restraint_wt  =  50.0               )

save_PDB_to_file(system,os.path.join(PEPDICE_EXAMPLES, 'data/02_geo_opt.pdb'))





insert_fragment (molecule = system,
                            fragment   = fragment,
                            sidechain  = True)
pprint(fragment)

system.energy(log = True)



save_PDB_to_file(system,os.path.join(PEPDICE_EXAMPLES, 'data/03_fragment.pdb'))


'''
i = 1
monte_carlo(molecule           = system      ,
            #random             = random      ,
            temperature        = 10000          ,
            Kb                 = 1           ,
            angle_range        = 5           ,
            nSteps             = 5000        ,
            fragment_rate      = 1.0         , #between 0  and 1
            fragment_sidechain = True        ,
            PhiPsi_rate        = 0.0         ,
            trajectory         = TRAJECTORY+str(i),
            pn                 = 1                       )
'''
