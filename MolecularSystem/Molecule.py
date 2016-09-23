#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  Atom.py
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
import os
from Atom import Atom
from Coordinates import Coordinates
from Energy      import Energy
from Residue import  Residue
from Bio.PDB import PDBParser
import subprocess
'''
#!/bin/bash

# . Bash environment variables and paths to be added to a user's ".bash_profile" file.
# . Some of these values may need modifying (e.g. PEPDICE_SCRATCH and PYTHONPATH).

# . The root of the program.
PEPDICE_ROOT=/home/farminf/Programas/PepDice ; export PEPDICE_ROOT

# . Package paths.
PEPDICE_BABEL=$PEPDICE_ROOT/Babel                     ; export PEPDICE_BABEL
PEPDICE_CORE=$PEPDICE_ROOT/Core                       ; export PEPDICE_CORE
PEPDICE_MOLECULE=$PEPDICE_ROOT/MolecularSystem        ; export PEPDICE_MOLECULE
#PEPDICE_MOLECULESCRIPTS=$PEPDICE_ROOT/pMoleculeScripts-1.9.0 ; export PEPDICE_PMOLECULESCRIPTS

# . Additional paths.
PEPDICE_PARAMETERS=$PEPDICE_ROOT/Parameters                                   ; export PEPDICE_PARAMETERS
PEPDICE_SCRATCH=$PEPDICE_ROOT/scratch                                         ; export PEPDICE_SCRATCH
#PEPDICE_STYLE=$PEPDICE_PARAMETERS/ccsStyleSheets/defaultStyle.css ; export PEPDICE_STYLE

# . The python path.
PYTHONPATH=:$PEPDICE_ROOT/Babel:$PEPDICE_ROOT/Core:$PEPDICE_ROOT/MolecularSystem ; export PYTHONPATH
'''


def pdb_corrections (pdbin = None, pdbout = None):
    """ Function doc """
    pdb  = open(pdbin, 'r') 
    pdb  = pdb.readlines()
    print pdb
    text = []
    #for line in pdb:
    #    print line


#-------------------------------------------------------------------------------
PEPDICE = os.environ.get('PEPDICE')
PEPDICE_EXAMPLES = os.path.join(PEPDICE, 'Examples')
PEPDICE_PARAMETER= os.path.join(PEPDICE, 'Parameters')
#-------------------------------------------------------------------------------
from AATorsions import load_torsion_from_file

class Molecule(Atom   ,
               Residue,
               Coordinates,
               Energy
               ):
    """ class to store info about a molecule"""

    def __init__(self, id=0, name = 'protein'):

        self.id       = id
        self.name     = name
        self.residues = []

        self.top     = None
        self.psf     = None
        self.param   = None
        self.ff_type = None


        self.concections = {
                            'bonds': [],
                            'angles': [],
                            'dihedral': [],
                            'improper': [],
                           }

        self.torsions       = None
        self.FIX_atoms_CHI  = None
        # Parameters
        self.fixed_residues = []
        self.fragments      = []
        
        
        
        # MC important atributes
        self.actual_energy   = None
        self.previous_energy = None 
        
        
        
        
        

        #self.pn
        self.bond      = 1.0
        self.angle     = 1.0
        self.dihed     = 1.0
        self.imprp     = 1.0
        self.elect     = 1.0
        self.vdw       = 1.0
        self.boundary  = 1.0
        self.esurf     = 1.0
        self.egb       = 1.0
        self.AB        = 0.0

    def build_peptide_from_sequence_AMBER (self, sequence = None, force_field = 'ff03ua', overwrite   = True  , NCTER = False):
        """ 
        Function doc 
        source leaprc.ff03ua
        foo = sequence { ACE ALA NME }
        saveamberparm foo foo.top foo.crd
        """
        
        aa_dic = { 
                 'A' : 'ALA',
                 'R' : 'ARG',
                 'N' : 'ASN',
                 'D' : 'ASP',
                 'C' : 'CYS',
                 'E' : 'GLU',
                 'Q' : 'GLN',
                 'G' : 'GLY',
                 'H' : 'HIS',
                 'I' : 'ILE',
                 'L' : 'LEU',
                 'K' : 'LYS',
                 'M' : 'MET',
                 'F' : 'PHE',
                 'P' : 'PRO',
                 'S' : 'SER',
                 'T' : 'THR',
                 'W' : 'TRP',
                 'Y' : 'TYR',
                 'V' : 'VAL'
                 }
        
        if overwrite:
            text  = 'addpath '+PEPDICE+'/Parameters/amber/labio.amber \n'
            text  += 'source leaprc.' + force_field + ' \n'
            
            if NCTER:
                text += 'foo = sequence {N'
            else:
                text += 'foo = sequence {'

            n = 1
            n2 = 1
            for aa in sequence:
                
                if n == len(sequence):
                    if NCTER:
                        text += 'C'
                    else:
                        pass
                if n2 >= 10:
                    text += aa_dic[aa] + ' \n'
                    n2 = 1

                else:
                    text += aa_dic[aa] + ' '
                    n2 += 1
                n += 1
                
                
                
            text += '} \n'
            text += 'saveamberparm foo '+ self.name +'.top '+ self.name +'.crd \n'
            text += 'savepdb foo ' + self.name +'.pdb \n'
            text += 'quit'
            leaprc = open('leaprc', 'w')
            leaprc.write(text)
            leaprc.close()
            
        
        
        
        
        #os.system('tleap -f leaprc')
        subprocess.call(['tleap', '-f', 'leaprc'])
        pdb_corrections (pdbin = self.name+'.pdb', pdbout = None)
        #subprocess.call(['babel', '-ipdb', self.name +'.pdb' ,'-opdb', self.name +'.pdb'])
        
        self.load_PDB_to_system      (filename  = self.name + '.pdb')  
        self.import_AMBER_parameters (top       = self.name + '.top' ,
                                      torsions  = os.path.join(PEPDICE_PARAMETER, 'amber/AMBER_rotamers.dat') )

    def build_peptide_from_sequence (self,
                                     sequence    = 'AAAAAAAAA',
                                     _type       = 'amber'    ,
                                     force_field = 'ff03ua'    ,
                                     overwrite   = True       ,
                                     ):
        """ :P """
        if _type == 'amber':
            self.build_peptide_from_sequence_AMBER(sequence    = sequence   , 
                                                   force_field = force_field, 
                                                   overwrite   = overwrite  )





    def load_PDB_to_system(self, filename = None):
        parser    = PDBParser(QUIET=True)
        structure = parser.get_structure('X', filename)
        self.residues = []

        for model in structure:

            c = 1
            for chain in model:

                self.id   = 1
                #self.name = "protein"

                n = 1
                r = 1

                for pdb_residue in chain:
                    residue = Residue(id=r,  name=pdb_residue.resname)
                    for pdb_atom in pdb_residue:

                        atom = Atom(id=n,
                                    name=pdb_atom.name,
                                    pos=pdb_atom.coord)
                        n += 1

                        residue.atoms.append(atom)
                    self.residues.append(residue)
                    r += 1


    def import_CHARMM_parameters (self, psf = None, param = None, torsions = None):
        """ Function doc """
        self.psf          = psf
        self.param        = param
        self.ff_type      = 'charmm'
        self.torsions     = load_torsion_from_file (torsions)

        sef.FIX_atoms_CHI ={
                                'CHI1' : ['OT1','OT2','CA','N','H', 'HT1', 'HT2', 'HT3', 'HN','C', 'O','HA'  ],

                                'CHI2' : ['OT1','OT2','CA','N','H', 'HT1', 'HT2', 'HT3', 'HN','C', 'O','HA',
                                          'HB','HB1','HB2','CG2', 'HG21', 'HG22', 'HG23'         ],

                                'CHI3' : ['OT1','OT2','CA','N','H', 'HT1', 'HT2', 'HT3', 'HN','C', 'O','HA',
                                          'HB','HB1','HB2', 'CB', 'HG1', 'HG2'                   ],

                                'CHI4' : ['OT1','OT2','CA','N','H', 'HT1', 'HT2', 'HT3', 'HN','C', 'O','HA',
                                          'HB','HB1','HB2', 'CB', 'HG1', 'HG2', 'CG', 'HD1','HD2'],

                                'CHI5' : ['OT1','OT2','CA','N','H', 'HT1', 'HT2', 'HT3', 'HN','C', 'O','HA',
                                          'HB','HB1','HB2', 'CB', 'HG1', 'HG2', 'CG', 'HD1','HD2',
                                          'HE' ,'HE1', 'HE2', 'CD'                               ],
                               }


    def import_AMBER_parameters (self, top = None, torsions = None):
        """ Function doc """
        self.top          = top
        self.ff_type      = 'amber'
        self.torsions     = load_torsion_from_file (torsions)
        self.FIX_atoms_CHI ={
                                'CHI1' : ['OT1','OT2','CA','N','H', 'H1', 'H2', 'H3', 'HN','C', 'O','HA'  ],

                                'CHI2' : ['OT1','OT2','CA','N','H', 'H1', 'H2', 'H3', 'HN','C', 'O','HA',
                                          'CB','HB', 'HB2', 'HB3', 'CG2', 'HG21', 'HG22', 'HG23',
                                         ],

                                         # 'HB','HB1','HB2','HB3', 'CG2', 'HG21', 'HG22','HG23','HG3'      ],

                                'CHI3' : ['OT1','OT2','CA'  ,'N'  ,'H'  ,'H1'  ,'H2'  ,'H3'  ,'HN','C','O','HA',
                                          'CB' ,'HB' , 'HB2','HB3','CG2','HG21','HG22','HG23',
                                          'CG' ,'HG2', 'HG3',
                                          ],

                                'CHI4' : ['OT1','OT2','CA','N','H', 'H1', 'H2', 'H3', 'HN','C', 'O','HA' ,
                                          'CB' ,'HB' , 'HB2','HB3','CG2','HG21','HG22','HG23',
                                          'CG' ,'HD1','HG2','HG3',
                                          'CD', 'HD1','HD2','HD3'
                                          ],

                                'CHI5' : ['OT1','OT2','CA','N','H', 'H1', 'H2', 'H3', 'HN','C', 'O','HA' ,
                                          'CB' ,'HB' , 'HB2','HB3','CG2','HG21','HG22','HG23',
                                          'CG' ,'HG1','HG2','HG3' ,
                                          'CD' , 'HD1','HD2','HD3',
                                          'HE' ,'HE1','HE2','HE3'
                                          ],
                               }


    def import_Calpha_model_parameters (self, top = None, torsions = None):
        """ Function doc """
        #self.top          = top
        self.ff_type      = 'Calpha_model'

        #----------------------------------------------------------------#
        #           Atribuicao dos valores de Hydropathy Index           #
        #----------------------------------------------------------------#
                                                                         #
        '''                                                              #
        To test hydrophobicity within sub-conformations we have used     #
        Kyte’s and Doolittle’s ‘Hydropathy Index’ [10] as a model.       #
        This works by assigning each amino acid a value to represents    #
        its hydrophobic (H) or hydrophilic (P) properties. The larger    #
        the value is the more hydrophobic the amino acid is(see Table 1).#
                                                                         #
                        Table 1.                                         #
        Kyte’s and Doolittle’s Hydropathy Index [10]                     #
        Amino Acid   Three-Letter  Hydropathy Index                      #
        Glycine          GLY           -0.4                              #
        Alanine          ALA            1.8                              #
        Proline          PRO            1.6                              #
        Valine           VAL            4.2                              #
        Leucine          LEU            3.8                              #
        Isoleucine       ILE            4.5                              #
        Methionine       MET            1.9                              #
        Phenylalanine    PHE            2.8                              #
        Tyrosine         TYR           -1.3                              #
        Tryptophan       TRP           -0.9                              #
        Serine           SER           -0.8                              #
        Threonine        THR           -0.7                              #
        Cysteine         CYS            2.5                              #
        Asparagine       ASN           -3.5                              #
        Glutamine        GLN           -3.5                              #
        Lysine           LYS           -3.9                              #
        Histidine        HIS           -3.2                              #
        Arginine         ARG           -4.5                              #
        Aspartate        ASP           -3.5                              #
        Glutamate        GLU           -3.5                              #
        '''                                                              #
                                                                         #
                                                                         #
        hydropathy_index = {                                             #
                            'GLY': -0.4,                                 #
                            'ALA':  1.8,                                 #
                            'PRO':  1.6,                                 #
                            'VAL':  4.2,                                 #
                            'LEU':  3.8,                                 #
                            'ILE':  4.5,                                 #
                            'MET':  1.9,                                 #
                            'PHE':  2.8,                                 #
                            'TYR': -1.3,                                 #
                            'TRP': -0.9,                                 #
                            'SER': -0.8,                                 #
                            'THR': -0.7,                                 #
                            'CYS':  2.5,                                 #
                            'ASN': -3.5,                                 #
                            'GLN': -3.5,                                 #
                            'LYS': -3.9,                                 #
                            'HIS': -3.2,                                 #
                            'ARG': -4.5,                                 #
                            'ASP': -3.5,                                 #
                            'GLU': -3.5                                  #
                            }                                            #
                                                                         #
        for residue_i in self.residues:                                  #
            for atom_i in residue_i.atoms:                               #
                if atom_i.name == 'CA':                                  #
                    name = residue_i.name                                #
                    if name in ['HSD', 'HSE', 'HDP', 'HIE', 'HID']:      #
                        name = 'HIS'                                     #
                    atom_i.epsilon = hydropathy_index[name]              #


    def import_GMX_parameters (self, tpr = None, torsions = None):
        """ Function doc """
        self.tpr          = tpr
        self.ff_type      = 'gmx'
        #self.torsions     = load_torsion_from_file (torsions)
        self.FIX_atoms_CHI ={
                                'CHI1' : ['OT1','OT2','CA','N','H', 'H1', 'H2', 'H3', 'HN','C', 'O','HA'  ],

                                'CHI2' : ['OT1','OT2','CA','N','H', 'H1', 'H2', 'H3', 'HN','C', 'O','HA',
                                          'CB','HB', 'HB2', 'HB3', 'CG2', 'HG21', 'HG22', 'HG23',
                                         ],

                                         # 'HB','HB1','HB2','HB3', 'CG2', 'HG21', 'HG22','HG23','HG3'      ],

                                'CHI3' : ['OT1','OT2','CA'  ,'N'  ,'H'  ,'H1'  ,'H2'  ,'H3'  ,'HN','C','O','HA',
                                          'CB' ,'HB' , 'HB2','HB3','CG2','HG21','HG22','HG23',
                                          'CG' ,'HG2', 'HG3',
                                          ],

                                'CHI4' : ['OT1','OT2','CA','N','H', 'H1', 'H2', 'H3', 'HN','C', 'O','HA' ,
                                          'CB' ,'HB' , 'HB2','HB3','CG2','HG21','HG22','HG23',
                                          'CG' ,'HD1','HG2','HG3',
                                          'CD', 'HD1','HD2','HD3'
                                          ],

                                'CHI5' : ['OT1','OT2','CA','N','H', 'H1', 'H2', 'H3', 'HN','C', 'O','HA' ,
                                         'CB' ,'HB' , 'HB2','HB3','CG2','HG21','HG22','HG23',
                                          'CG' ,'HG1','HG2','HG3',
                                          'CD', 'HD1','HD2','HD3',
                                          'HE' ,'HE1','HE2','HE3'
                                          ],
                               }


    def import_fixed_from_string(self, fixed=None):
        """ Function doc """

        #-------------------importing from porter------------------------------
        n = 0
        for aa in fixed:
            if aa != '0':
                self.fixed_residues.append(n)
            n += 1
        #----------------------------------------------------------------------

    def PrintStatus (self, parameters = True, coordinates = True ):
        """ Function doc """
        print 'Bond_Stretching  = ', self.Bond_Stretching
        print 'Angle_Bending    = ', self.Angle_Bending
        print 'Improper_Torsion = ', self.Improper_Torsion
        print 'Torsional_Angle  = ', self.Torsional_Angle
        print 'Van_der_Waals    = ', self.Van_der_Waals
        print 'Charge_Charge    = ', self.Charge_Charge
        print 'Total_Energy     = ', self.Total_Energy

        print 'Fixed_residues   = ', self.fixed_residues

        for res in self.residues:
            for atom in res.atoms:
                print atom.name, atom.pos


    def import_restraints_from_porter (self, filein= None):
        """ Function doc """
        text = open(filein, 'r')
        text =  text.readlines()

        for line in text:
            line2 =  line.split()
            if line2[0] == 'rest':
                fixed_residues = line2[2]

        for resi in fixed_residues:
            try:
                self.fixed_residues.append(int(resi))
            except:
                pass
        print len(self.fixed_residues), self.fixed_residues

