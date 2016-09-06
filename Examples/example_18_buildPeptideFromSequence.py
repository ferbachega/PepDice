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
#from Molecule    import Molecule                     #
from Geometry    import *                            #
from MonteCarlo  import monte_carlo, insert_fragment                  #
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
from Fragments import import_fragments_from_pdb


#eita

#-------------------------------------------------------------------------------
PEPDICE = os.environ.get('PEPDICE')
PEPDICE_EXAMPLES = os.path.join(PEPDICE, 'Examples')
PEPDICE_PARAMETER= os.path.join(PEPDICE, 'Parameters')
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
#system = Molecule() 
#system.load_PDB_to_system      (filename = os.path.join(PEPDICE_EXAMPLES , 'data/alpha/1GAB/1gab_ff03ua_AMBER_extended.pdb'   )   )   
#system.import_AMBER_parameters(top       = os.path.join(PEPDICE_EXAMPLES , 'data/alpha/1GAB/1gab_ff03ua_AMBER_extended.prmtop')   ,   
#                                torsions = os.path.join(PEPDICE_PARAMETER, 'amber/AMBER_rotamers.dat') )   
#TRAJECTORY  = '/home/farminf/Programas/PepDice/Examples/1GAB/amber/1gab_amber_side_chain_rand.xyz' /home/farminf/Programas/pepdice/Examples/data/alpha/1GAB/1gab_ff03ua_AMBER_extended.prmtop
#-------------------------------------------------------------------------------


system = Molecule() 
system.name =  'teste_novos_decoys'
system.build_peptide_from_sequence ( sequence    = 'GSPEVQIAILTEQINNLNEHLRVHKKDHHSRRGLLKMVGKRRRLLAYLRNKDVARYREIVEKLGL',
                                     _type       = 'amber'    ,
                                     force_field = 'ff03ua'   ,
                                     overwrite   = True       ,
                                     )
#---------------------------------------------------------------------------------


#---------------------------------------------------------------------------------
#import_SS_from_string (molecule = system, ss ='CCHHHHHHHHHHHHHHHHHHHHHCCCCCHHHHHHHHHHHHHHHHHHHHHHHCHHHHHHHHHHHCC')
#save_PDB_to_file(system, 'teste_novos_decoys_SS.pdb')
#---------------------------------------------------------------------------------



#---------------------------------------------------------------------------------
pdb     = 'native.pdb'      
system2 = Molecule()                                                                                                                 #
system2.load_PDB_to_system(filename = pdb)

fragment = import_fragments_from_pdb (molecule  = system2, 
                                      residues  = range(0,64), 
                                      mainchain = True, 
                                      sidechain = False)
#---------------------------------------------------------------------------------





#---------------------------------------------------------------------------------
#minimize(molecule = system,
#               imin  = 1          ,
#               maxcyc= 1000        ,
#               ncyc  = 100        ,
#               cut   = 10         ,
#               rgbmax= 999        ,
#               igb   = 1          ,
#               ntb   = 0          ,
#               ntpr  = 100        ,
#               ntr   = 0          )
#               #restraintmask = ':1-50 & @CA,N,C,O=', 
#               #restraint_wt  =  50.0 
#
#---------------------------------------------------------------------------------

insert_fragment (molecule = system,
                            fragment   = fragment,
                            sidechain  = True)
system.energy(log = True)

save_PDB_to_file(system,os.path.join(PEPDICE_EXAMPLES, 'teste_novos_decoys_fragment.pdb'))


















print system.energy( log = True)

