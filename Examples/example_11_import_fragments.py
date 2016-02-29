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
from SideChainRefine import optimize_side_chain
#----------------------------------------------------#

#-------------------------------------------------------------------------------
PEPDICE = os.environ.get('PEPDICE')
PEPDICE_EXAMPLES = os.path.join(PEPDICE, 'Examples')
PEPDICE_PARAMETER= os.path.join(PEPDICE, 'Parameters')
#-------------------------------------------------------------------------------
#/home/farminf/Programas/PepDice/Examples/outputs/1gab_amber_example04_extended.pdb

#---------------------------------------------------------------------------------------------------------
system = Molecule() 
system.load_PDB_to_system      (filename = os.path.join(PEPDICE_EXAMPLES , 'data/1gab_amber03ua.pdb'   )   )  
system.import_AMBER_parameters(top       = os.path.join(PEPDICE_EXAMPLES , 'data/1gab_amber03ua.prmtop')   ,   
                                torsions = os.path.join(PEPDICE_PARAMETER, 'amber/AMBER_rotamers.dat') )   

TRAJECTORY  = os.path.join(PEPDICE_EXAMPLES , 'outputs/1gab_amber_example10_monte_carlo.xyz')
#---------------------------------------------------------------------------------------------------------
try:
    os.remove(TRAJECTORY)
except:
    pass
    
    
    

#system.energy(log = True)

def import_fragments_from_pdb (molecule = None, residues = [], mainchain =True, sidechain = True):
    """ Function doc """
    fragment = {}
    
    for i in residues:
        if mainchain:
            fragment[i] = {}
            fragment[i]['PHI']   = computePhiPsi (molecule=molecule, resi=i, bond='PHI')
            fragment[i]['PSI']   = computePhiPsi (molecule=molecule, resi=i, bond='PSI')
            fragment[i]['OMEGA'] = computePhiPsi (molecule=molecule, resi=i, bond='OMEGA')
            fragment[i]['NAME']  = molecule.residues[i].name
        if sidechain:
            for chi in ["CHI1","CHI2","CHI3","CHI4","CHI5"]:
                fragment[i][chi] = computeCHI (molecule= molecule, resi=i, bond=chi)
                #print i, system.residues[i].name, chi,computeCHI (molecule=system, resi=i, bond=chi)
    return fragment 





molecule_temp = Molecule()
molecule_temp.load_PDB_to_system(filename = os.path.join(PEPDICE_EXAMPLES ,'data/1gab_files/1gab_0001.pdb'))
molecule_temp.import_AMBER_parameters(top   = os.path.join(PEPDICE_EXAMPLES , 'data/1gab_amber03ua.prmtop')   ,   
                                   torsions = os.path.join(PEPDICE_PARAMETER, 'amber/AMBER_rotamers.dat') ) 

print os.path.join(PEPDICE_EXAMPLES ,'data/1gab_files/1gab_0001.pdb')
print len(molecule_temp.residues)


pdbs          = [
                 os.path.join(PEPDICE_EXAMPLES ,'data/1gab_files/1gab_0001.pdb'),
                 os.path.join(PEPDICE_EXAMPLES ,'data/1gab_files/1gab_0002.pdb'),
                 os.path.join(PEPDICE_EXAMPLES ,'data/1gab_files/1gab_0003.pdb'),
                 os.path.join(PEPDICE_EXAMPLES ,'data/1gab_files/1gab_0004.pdb'),
                 os.path.join(PEPDICE_EXAMPLES ,'data/1gab_files/1gab_0005.pdb'),
                 os.path.join(PEPDICE_EXAMPLES ,'data/1gab_files/1gab_0006.pdb'),
                 os.path.join(PEPDICE_EXAMPLES ,'data/1gab_files/1gab_0007.pdb'),
                 os.path.join(PEPDICE_EXAMPLES ,'data/1gab_files/1gab_0008.pdb'),
                 os.path.join(PEPDICE_EXAMPLES ,'data/1gab_files/1gab_0009.pdb'),
                 os.path.join(PEPDICE_EXAMPLES ,'data/1gab_files/1gab_0010.pdb'),
                 os.path.join(PEPDICE_EXAMPLES ,'data/1gab_files/1gab_0011.pdb'),
                 os.path.join(PEPDICE_EXAMPLES ,'data/1gab_files/1gab_0012.pdb'),
                 os.path.join(PEPDICE_EXAMPLES ,'data/1gab_files/1gab_0013.pdb'),
                 os.path.join(PEPDICE_EXAMPLES ,'data/1gab_files/1gab_0014.pdb'),
                 os.path.join(PEPDICE_EXAMPLES ,'data/1gab_files/1gab_0015.pdb'),
                 os.path.join(PEPDICE_EXAMPLES ,'data/1gab_files/1gab_0016.pdb'),
                 os.path.join(PEPDICE_EXAMPLES ,'data/1gab_files/1gab_0017.pdb'),
                 os.path.join(PEPDICE_EXAMPLES ,'data/1gab_files/1gab_0018.pdb'),
                 os.path.join(PEPDICE_EXAMPLES ,'data/1gab_files/1gab_0019.pdb'),
                 os.path.join(PEPDICE_EXAMPLES ,'data/1gab_files/1gab_0020.pdb')]

#
pdbs          = [os.path.join(PEPDICE_EXAMPLES ,'data/1gab_files/2bjb.pdb')]
#
#

system.fragments = []
frag_size        = 7

for pdb in pdbs:
    molecule_temp.load_PDB_to_system(filename = pdb)  
    for i in range(0,100):
        resi     = random.randint(0, len(system.residues)-frag_size)
        #print resi
        fragment = import_fragments_from_pdb (molecule = molecule_temp, 
                                              residues = range(resi, resi+frag_size), 
                                             sidechain = True)
        #pprint (fragment)
        system.fragments.append(fragment)
pprint(system.fragments)
#print len(system.fragments)





#for i in range(0,100):
#    resi     = random.randint(0, len(system.residues)-frag_size)
#    fragment = import_fragments_from_pdb (molecule = system, 
#                                          residues = range(resi,resi+frag_size), 
#                                         sidechain = True)
#    #pprint (fragment)
#    system.fragments.append(fragment)
##pprint(system.fragments)
#
#
#system.load_PDB_to_system (filename = os.path.join(PEPDICE_EXAMPLES , 'outputs/1gab_amber_example04_ff03ua_extended.pdb'))
#monte_carlo(molecule   = system      ,
#           temperature = 300         ,
#           Kb          = 0.0083144621,
#           angle_range = 45          ,
#           nSteps      = 1000        ,
#           trajectory  = TRAJECTORY)
#
