#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  Fragments.py
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
from random import randint                                                                          #
                                                                                                    #
from Energy import save_PDB_to_file                                                                 #
                                                                                                    #
#---------------------------------------------------------------------------------------------------#





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


def build_fragment_library_from_pdbs (
                                     molecule             = None ,
                                     frag_size            = 3    ,
                                     number_of_fragments  = 100  ,
                                     pdblist              = []   ,
                                     ):
    """ Function doc """
    for pdb in pdblist:
        molecule.load_PDB_to_system(filename = pdb)  
        for i in range(0,number_of_fragments):
            resi     = random.randint(0, len(molecule.residues)-frag_size)
            fragment = import_fragments_from_pdb (molecule = molecule, 
                                                  residues = range(resi, resi+frag_size), 
                                                 sidechain = True)
            molecule.fragments.append(fragment)
    return molecule



