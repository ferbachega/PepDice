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
import pprint
import random


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
                residue   = molecule.residues[i]
                if chi in molecule.torsions[residue.name]:
                    #print chi
                    try:
                        fragment[i][chi] = computeCHI (molecule= molecule, resi=i, bond=chi)
                    except:
                        print i , 'and', chi ,  ' failed'
                
                #print i, system.residues[i].name, chi,computeCHI (molecule=system, resi=i, bond=chi)
    return fragment 


def build_fragment_library_from_pdbs (
                                     molecule             = None ,
                                     frag_size            = 3    ,
                                     number_of_fragments  = 100  ,
                                     pdblist              = []   ,
                                     ):
    """ Function doc """
    
    
    
    fragments =  []
    '''
    fragments =  [
                  [     -> lista  indica a posicao  na sequencia alvo
                   
                   {},  -> dada uma posicao, exista N possiveis fragmentos {}
                   {}, ...
                  
                  ],
                  
                  [],
                  
                  [],
                 ]
    '''

    #print 1, len(molecule.residues)
    # lista  (posicoes  =  index do residuo), 
    # cade elemento eh uma lista com N fragmentos 
    # possiveis.
    for resi in range(len(molecule.residues)-frag_size): 
        fragments.append([])
    
    
    

    for pdb in pdblist:
        #para os pdbs na lista, importar fragmentos
        molecule.load_PDB_to_system(filename = pdb)  
        #print 2, len(molecule.residues), pdb
        
        for i in range(0, number_of_fragments):
            
            # sortear uma posicao na sequencia target
            resi     = random.randint(0, len(molecule.residues)-frag_size -1)
            
            # importar o fragment  resi+frag_size - normalmente entre 3-9
            fragment = import_fragments_from_pdb (molecule = molecule, 
                                                  residues = range(resi, resi+frag_size), 
                                                 sidechain = True)
            
            # adicionar o fragment na lista (resi = posicao do alinhamento) 
            fragments[resi].append(fragment)
            molecule.fragments = fragments
    
    print len(fragments)
    print len(fragments[0])
    n = 0 
    for resi in fragments:
        k = 0
        for frag in resi: 
            print 'Position: ',n , 'fragment index: ',k, 'Number of fragments : ',len(resi), 'fragment size : ', len(frag)
            k += 1
        n += 1
    
    return molecule



