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
#from MonteCarlo  import monte_carlo,  monte_carlo_dic, MC_replica_exchange, run_MC_replica_exchange #
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


'''

class Fragment:
    """ Class doc """
    
    def __init__ (self):
        """ Class initialiser """
        self.phi   = None
        self.psi   = None
        self.omega =  None
        
        pass

'''



def import_fragments_from_pdb (molecule = None, residues = [], mainchain =True, sidechain = True):
    """ Function doc """
    fragment = {}
    
    for i in residues:
        
        if mainchain:
            fragment[i] = {}
            fragment[i]['PHI']       = computePhiPsi (molecule=molecule, resi=i, bond='PHI')
            fragment[i]['PSI']       = computePhiPsi (molecule=molecule, resi=i, bond='PSI')
            fragment[i]['OMEGA']     = computePhiPsi (molecule=molecule, resi=i, bond='OMEGA')
            fragment[i]['NAME']      = molecule.residues[i].name
            #fragment[i]['position']  = None
            
            
        if sidechain:
        
            for chi in ["CHI1","CHI2","CHI3","CHI4","CHI5"]:
                residue   = molecule.residues[i]
                if chi in molecule.torsions[residue.name]:
                    #print i, residue.name, chi,  molecule.torsions[residue.name]
                    try:
                        fragment[i][chi] = computeCHI (molecule= molecule, resi=i, bond=chi)
                        #print i, residue.name, chi,  fragment[i][chi] 
                    except:
                        fragment[i][chi] = None
                        #print i, residue.name, chi,  fragment[i][chi] 
                
                #print i, system.residues[i].name, chi,computeCHI (molecule=system, resi=i, bond=chi)
    return fragment 


def build_fragment_library_from_pdbs (
                                     molecule             = None ,
                                     frag_size            = 3    ,
                                     number_of_fragments  = 100  ,
                                     pdblist              = []   ,
                                     log                  = True ,
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

    
    

    fileout = open ('fragment_list.log', 'w')
    text    = 'seq:   '
    
    #print len(molecule.residues*'X')
    #print 1, len(molecule.residues)
    # lista  (posicoes  =  index do residuo), 
    # cade elemento eh uma lista com N fragmentos 
    # possiveis.
    
    for resi in range(len(molecule.residues)-frag_size): 
        fragments.append([])
    
    for resi in range(len(molecule.residues)): 
        #fragments.append([])
        #text += ('X')
        text += molecule.AminoAcid_dic[molecule.residues[resi].name]
    
    text+='\n'
    
    

    for pdb in pdblist:
        #para os pdbs na lista, importar fragmentos
        molecule.load_PDB_to_system(filename = pdb)  
        #print 2, len(molecule.residues), pdb
        
        for i in range(0, number_of_fragments):
            text += '%-7s' %(i)
            # sortear uma posicao na sequencia target
            resi     = random.randint(0, len(molecule.residues)-frag_size -1)
            
            # importar o fragment  resi+frag1_size - normalmente entre 3-9
            fragment = import_fragments_from_pdb (molecule = molecule, 
                                                  residues = range(resi, resi+frag_size), 
                                                 sidechain = True)
            
            
            text_fragment = '  '
            
            position  = ''
            for k in range(len(molecule.residues)):
                if k in fragment:
                    text    += 'X' 
                    position+= 'X'
                    
                    if fragment[k]['PHI'] == None:
                        print     'failed (phi):',k, fragment[k]['PHI']
                        fragment[k]['PHI'] = 0.0

                    if fragment[k]['PSI'] == None:
                        print     'failed (psi):',k, fragment[k]['PSI']
                    
                        fragment[k]['PSI'] = 0.0
                    #text_fragment += 'k %5d phi:%8.5f   psi:%8.5f  omega:' %(k, fragment[k]['PHI'],fragment[k]['PSI'])
                    text_fragment += 'position: %-5d  phi:%14.5f psi:%14.5f omega:%14.5f  |  ' %(k, fragment[k]['PHI'],fragment[k]['PSI'],fragment[k]['OMEGA'] )

                else:
                    text    += '-'
                    position+= '-'
           
            
            #fragment['position'] = position
            #print fragment['position']
            

            text += text_fragment
            text += '\n'
            
            # adicionar o fragment na lista (resi = posicao do alinhamento) 
            fragments[resi].append(fragment)
            molecule.fragments = fragments
        
        text += '\n'
        fileout.write(text)

    
    
    #print len(fragments)
    #print len(fragments[0])
    n = 0 
    
    for resi in fragments:
        k = 0
        for frag in resi: 
            print 'Position: ',n , 'fragment index: ',k, 'Number of fragments : ',len(resi), 'fragment size : ', len(frag)
            k += 1
        n += 1
    
    return molecule



