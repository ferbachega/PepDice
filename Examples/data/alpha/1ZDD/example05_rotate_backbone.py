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




#-------------------------------------------------------------------------------
PEPDICE = os.environ.get('PEPDICE')
PEPDICE_EXAMPLES = os.path.join(PEPDICE, 'Examples')
PEPDICE_PARAMETER= os.path.join(PEPDICE, 'Parameters')
#-------------------------------------------------------------------------------


#---------------------------------------------------------------------------------------------------------
system = Molecule() 
system.load_PDB_to_system      (filename = os.path.join(PEPDICE_EXAMPLES , 'data/1gab_amber.pdb'   )   )   
system.import_AMBER_parameters(top       = os.path.join(PEPDICE_EXAMPLES , 'data/1gab_amber.prmtop')   ,   
                                torsions = os.path.join(PEPDICE_PARAMETER, 'amber/AMBER_rotamers.dat') )   

TRAJECTORY  = os.path.join(PEPDICE_EXAMPLES , 'outputs/1gab_amber_example04_refold.xyz')
#---------------------------------------------------------------------------------------------------------


def test_rotetaPhiPsiOMEGA (system = None, theta =  0.017453292519943295,  computeTorsion = False):
    """ Function doc """
    backup  = system 
    
    for j in range(0,10):
        system = backup
        for i in range (0,len(system.residues)):

            try:
                rotate_backbone(molecule=system, resi=i, bond='PHI'  , theta= theta*(randint(j*-1,j)))
            except:
                print 'impossible to rotate PHI'
            
            try:
                rotate_backbone(molecule=system, resi=i, bond='PSI'  , theta=theta*(randint(j*-1,j)))
            except:
                print 'impossible to rotate PSI'

            try:
                rotate_backbone(molecule=system, resi=i, bond='OMEGA', theta=theta*(randint(j*-1,j)))
            except:
                print 'impossible to rotate OMEGA'
                
                
        save_PDB_to_file(system,  os.path.join(PEPDICE_EXAMPLES ,'outputs/1gab_amber_example05_rotetaBackbone_'+str(j)+'_rand.pdb'))
        #test_computeTorsions (system = system)

test_rotetaPhiPsiOMEGA (system = system)














