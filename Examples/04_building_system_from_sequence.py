#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  04_building_system_from_sequence.py
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
from MonteCarlo  import monte_carlo#,  monte_carlo_dic, MC_replica_exchange, run_MC_replica_exchange #
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
import random
random.seed(1234)


#-------------------------------------------------------------------------------
PEPDICE = os.environ.get('PEPDICE')
PEPDICE_EXAMPLES = os.path.join(PEPDICE, 'Examples')
PEPDICE_PARAMETER= os.path.join(PEPDICE, 'Parameters')
#-------------------------------------------------------------------------------
#/home/farminf/Programas/PepDice/Examples/outputs/1gab_amber_example04_extended.pdb

system = Molecule()

system.build_peptide_from_sequence (
                                     sequence    = 'AADD',
                                     _type       = 'amber'       ,
                                     force_field = 'ff03ua.labio',
                                     overwrite   = True          ,
                                     )
system.set_energy_model('FULL')

#system.energy(log = True)


hydropathic_table_AB = {
                            'ARG' : 'B', #-4.5,  
                            'LYS' : 'B', #-3.9,
                            'ASN' : 'B', #-3.5,
                            'ASP' : 'B', #-3.5,
                            'GLU' : 'B', #-3.5,
                            'GLN' : 'B', #-3.5,
                            'HIS' : 'B', #-3.2,
                            'HIE' : 'B', #-3.2,

                            'PRO' : 'B', #-1.6,
                            'TYR' : 'B', #-1.3,
                            'TRP' : 'B', #-0.9,
                            'SER' : 'B', #-0.8,
                            'THR' : 'B', #-0.7,
                            'GLY' : 'B', #-0.4,
                            'ALA' : 'A', # 1.8,
                            'MET' : 'A', # 1.9,
                            'CYS' : 'A', # 2.5,
                            'PHE' : 'A', # 2.8,
                            'LEU' : 'A', # 3.8,
                            'VAL' : 'A', # 4.2,
                            'ILE' : 'A', # 4.5,
                            }




name_i = system.residues[0].name
for atom in system.residues[0].atoms:
    if atom.name == 'CA':
        atom_i    = atom  
        atom_i.AB = 'A'
        #atom_i.sigma  = 3.8
        #mass_i = mass_table[name_i]
        
        

name_j = system.residues[1].name
for atom in system.residues[1].atoms:
    if atom.name == 'CA':
        atom_j = atom  
        atom_j.AB = 'A'
        #atom_j.sigma = 3.8

name_k = system.residues[2].name
for atom in system.residues[2].atoms:
    if atom.name == 'CA':
        atom_k = atom  
        atom_k.AB = 'B'
        #atom_j.sigma = 3.8

name_l = system.residues[3].name
for atom in system.residues[3].atoms:
    if atom.name == 'CA':
        atom_l = atom  
        atom_l.AB = 'B'
        #atom_j.sigma = 3.8



r = 3.5
for i in range (1, 600):
   r +=0.01
   
   print r, system.compute_atomi_atomj_AB_energy(atom_i, atom_j, r), system.compute_atomi_atomj_AB_energy(atom_l, atom_k, r), system.compute_atomi_atomj_AB_energy(atom_i, atom_k, r)
    



