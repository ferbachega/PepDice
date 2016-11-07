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
from Atom    import Atom
from ModelAB import ModelAB
from Coordinates import Coordinates
from Energy      import Energy
from Energy      import save_PDB_to_file
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

class Molecule(Atom       ,
               Residue    ,
               Coordinates,
               Energy     ,
               ModelAB    ,
               ):
    """ 
    class to store info about a molecule
    
    energy models:
        amber   - amber ff03ua energy components
        Calpha  - C alpha model  - like AB model
        Contact - Contact map only
        LPFSF   - LABIO protein folding scoring function
        
    """

    def __init__(self, id=0, name = 'protein'):
        
        self.energy_model  = 'amber'
        
        self.energy_models = {'amber'  ,
                              'Calpha' ,
                              'Contact',
                              'LSF'    ,
                              'FULL'   }
        
        self.id       = id
        self.name     = name
        self.residues = []

        self.top     = None
        self.psf     = None
        self.param   = None
        self.ff_type = None

        self.torsions       = None
        self.FIX_atoms_CHI  = None
        
        # Parameters
        self.fixed_residues = []
        self.fragments      = []
        self.cmap           = None 
        
        
        # MC important atributes
        self.actual_energy   = None
        self.previous_energy = None 
        
        
        # Energy components
        #self.pn
        
        
        self.energy_model_parameters = None
        
        self.R_contact = 0.0
        self.default_energy_model_setup ={'amber':{
                                                  #         AMBER
                                                  'ANGLE'      : 1.0       ,
                                                  'BOND'       : 1.0       ,
                                                  'DIHED'      : 1.0       ,
                                                  'EEL'        : 1.0       ,
                                                  'EELEC'      : 1.0       ,
                                                  'EGB'        : 1.0       ,
                                                  'EKtot'      : 1.0       ,
                                                  'EPtot'      : 1.0       ,
                                                  'ESURF'      : 1.0       ,
                                                  'Etot'       : 1.0       ,
                                                  'NB'         : 1.0       ,
                                                  'VDWAALS'    : 1.0       ,
                                                  
                                                  'cut'        : 999.0     ,
                                                  'igb'        : 1         ,
                                                  'saltcon'    : 0.2       ,
                                                  'gbsa'       : 1         ,
                                                  'rgbmax'     : 999.00000 ,
                                                  'surften'    : 0.010     ,
                                                  
                                                  'CONTACT'    : 0.0       ,
                                                  'R_contact'  : self.R_contact       ,
                                                  
                                                  'AB'  : 0.0       ,
                                                  'R_cutoff'   : 999.0     ,
                                                  
                                                  },

                                         'FULL':{
                                                  #         AMBER
                                                  'ANGLE'      : 1.0           ,
                                                  'BOND'       : 1.0           ,
                                                  'DIHED'      : 1.0           ,
                                                  'EEL'        : 1.0           ,
                                                  'EELEC'      : 1.0           ,
                                                  'EGB'        : 1.0           ,
                                                  'EKtot'      : 1.0           ,
                                                  'EPtot'      : 1.0           ,
                                                  'ESURF'      : 1.0           ,
                                                  'Etot'       : 1.0           ,
                                                  'NB'         : 1.0           ,
                                                  'VDWAALS'    : 1.0           ,
                                                                               
                                                  'cut'        : 999.0         ,
                                                  'igb'        : 1             ,
                                                  'saltcon'    : 0.2           ,
                                                  'gbsa'       : 1             ,
                                                  'rgbmax'     : 999.00000     ,
                                                  'surften'    : 0.010         ,
                                                                               
                                                  'CONTACT'    : 1.0           ,
                                                  'R_contact'  : self.R_contact,
                                                  
                                                  'AB'         : 1.0           ,
                                                  'R_cutoff'   : 999.0         ,
                                                  
                                                  },

                                          
                                          'Calpha' : {
                                                      'ANGLE'      : 0.0       ,
                                                      'BOND'       : 0.0       ,
                                                      'DIHED'      : 1.0       ,
                                                      'EEL'        : 0.0       ,
                                                      'EELEC'      : 0.0       ,
                                                      'EGB'        : 0.0       ,
                                                      'EKtot'      : 0.0       ,
                                                      'EPtot'      : 0.0       ,
                                                      'ESURF'      : 0.0       ,
                                                      'Etot'       : 0.0       ,
                                                      'NB'         : 0.0       ,
                                                      'VDWAALS'    : 1.0       ,
                                                      
                                                      'cut'        : 999.0     ,
                                                      'igb'        : 1         ,
                                                      'saltcon'    : 0.2       ,
                                                      'gbsa'       : 1         ,
                                                      'rgbmax'     : 999.00000 ,
                                                      'surften'    : 0.010     ,
                                                      
                                                      'CONTACT'    : 0.0       ,
                                                      'R_contact'  : self.R_contact       ,
                                                      
                                                      'AB'         : 1.0       ,
                                                      'R_cutoff'   : 999.0     ,
                                                  
                                                      },
                                          
                                          
                                          'Contact': {
                                                      'ANGLE'    : 0.0       ,
                                                      'BOND'     : 0.0       ,
                                                      'DIHED'    : 0.001     ,
                                                      'EEL'      : 0.0       ,
                                                      'EELEC'    : 0.0       ,
                                                      'EGB'      : 0.0       ,
                                                      'EKtot'    : 0.0       ,
                                                      'EPtot'    : 0.0       ,
                                                      'ESURF'    : 0.0       ,
                                                      'Etot'     : 0.0       ,
                                                      'NB'       : 0.0       ,
                                                      'VDWAALS'  : 0.001     ,
                                                      
                                                      'cut'        : 999.0     ,
                                                      'igb'        : 1         ,
                                                      'saltcon'    : 0.2       ,
                                                      'gbsa'       : 1         ,
                                                      'rgbmax'     : 999.00000 ,
                                                      'surften'    : 0.010     ,
                                                      
                                                      'CONTACT'    : 1.0       ,
                                                      'R_contact'  : self.R_contact       ,
                                                      
                                                      'AB'         : 0.0      ,
                                                      'R_cutoff'   : 0.0      ,
                                                  
                                                      },
                                          
                                          # energy = 1.15 -1.96E-5*energy_list['EEL'] -2.36E-5*energy_list['NB'] - 4.4E-4 *energy_list['DIHED'] + 1.85E-3*energy_list['VDWAALS'] - 7.5E-5*energy_list['EGB'] + 2.66E-5*energy_list['ESURF']

                                          'LSF'    :  {
                                                      
                                                      'info'       : ''' 
==============================================================================
Dep. Variable:                      y   R-squared:                       0.483
Model:                            OLS   Adj. R-squared:                  0.483
Method:                 Least Squares   F-statistic:                     1371.
Date:                Mon, 07 Nov 2016   Prob (F-statistic):               0.00
Time:                        19:32:29   Log-Likelihood:                -267.19
No. Observations:               13209   AIC:                             554.4
Df Residuals:                   13199   BIC:                             629.3
Df Model:                           9                                         
Covariance Type:            nonrobust                                         
==============================================================================
                 coef    std err          t      P>|t|      [95.0% Conf. Int.]
------------------------------------------------------------------------------
const          1.3119      0.017     75.956      0.000         1.278     1.346
SIZE          -0.0125      0.001    -13.212      0.000        -0.014    -0.011
CONTACT        0.0070      0.000     56.197      0.000         0.007     0.007
AB_ENERGY      0.3745      0.033     11.182      0.000         0.309     0.440
DIHED         -0.0003   5.22e-05     -6.478      0.000        -0.000    -0.000
EEL        -6.804e-05   7.37e-06     -9.236      0.000     -8.25e-05 -5.36e-05         
EELEC         -0.0002   1.04e-05    -15.772      0.000        -0.000    -0.000
EGB           -0.0001   1.01e-05    -13.817      0.000        -0.000    -0.000
ESURF          0.0317      0.001     60.724      0.000         0.031     0.033
VDWAALS        0.0479      0.003     15.507      0.000         0.042     0.054
==============================================================================
Omnibus:                      318.266   Durbin-Watson:                   1.694
Prob(Omnibus):                  0.000   Jarque-Bera (JB):              581.881
Skew:                          -0.189   Prob(JB):                    4.43e-127
Kurtosis:                       3.957   Cond. No.                     9.68e+04
==============================================================================
                                                                           ''' , 
                                                      
                                                      'CONSTANT'   :     1.3119,
                                                      'ANGLE'      :          0,
                                                      'BOND'       :          0,
                                                      'DIHED'      :    -0.0003,
                                                      'EEL'        : -6.804E-05,
                                                      'EELEC'      :    -0.0002,
                                                      'EGB'        :    -0.0001,
                                                      'EKtot'      :          0,
                                                      'EPtot'      :          0,
                                                      'ESURF'      :     0.0317,
                                                      'Etot'       :          0,
                                                      'NB'         :          0,
                                                      'VDWAALS'    :     0.0479,
                                                      'SIZE'       :    -0.0125,
                                                      
                                                      
                                                      
                                                      'cut'        : 999.0     ,
                                                      'igb'        : 1         ,
                                                      'saltcon'    : 0.2       ,
                                                      'gbsa'       : 1         ,
                                                      'rgbmax'     : 999.00000 ,
                                                      'surften'    : 0.010     ,
                                                      
                                                      'CONTACT'    :     0.0070,
                                                      'R_contact'  : self.R_contact,
                                                      
                                                      'AB'         :     0.3745,
                                                      'R_cutoff'   :      999.0,
                                                  
                                                      }
                                         }
                                         
        self.energy_model_parameters = self.default_energy_model_setup['amber']

                                         
        self.aa_dic = { 
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
    def set_energy_model (self, energy_model = 'amber'):
        """ Function doc """

        if energy_model in self.energy_models:
            self.energy_model  = energy_model
            self.energy_model_parameters = self.default_energy_model_setup[energy_model]

        else:
            print '\nEnergy model not found('+energy_model+'). Please select one of the options:'
            print self.energy_models
            print '\n'

    
    
    
    def build_peptide_from_sequence_AMBER (self, sequence = None, force_field = 'ff03ua', overwrite   = True  , NCTER = False):
        """ 
        Function doc 
        source leaprc.ff03ua
        foo = sequence { ACE ALA NME }
        saveamberparm foo foo.top foo.crd
        """
        

        
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
                    text += self.aa_dic[aa] + ' \n'
                    n2 = 1

                else:
                    text += self.aa_dic[aa] + ' '
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
                                     force_field = 'ff03ua'   ,
                                     overwrite   = True       ,
                                     ):
        """
         
        energy models:
            amber   - amber ff03ua energy components
            Calpha  - C alpha model  - like AB model
            Contact - Contact map only
            LPFSF   - LABIO protein folding scoring function
        
        """
        
        if _type == 'amber':
            self.build_peptide_from_sequence_AMBER(sequence    = sequence   , 
                                                   force_field = force_field, 
                                                   overwrite   = overwrite  )

        if _type == 'Calpha':
            pass
            
            
            input_sequence = ''
            for AA in sequence:
                if AA =='G':
                    input_sequence += 'G'
                else:
                    input_sequence += 'A'
                    
            self.build_peptide_from_sequence_AMBER(sequence    = input_sequence, 
                                                   force_field = force_field   ,
                                                   overwrite   = overwrite  )
            n = 0
            for res in self.residues:
                res.name = self.aa_dic[sequence[n]]
                n += 1
            save_PDB_to_file      (self, self.name + '.pdb')
            self.energy_model  = 'Calpha'
        
        
        if _type == 'contact':
            pass
            #self.build_peptide_from_sequence_AMBER(sequence    = sequence   , 
            #                                       force_field = force_field, 
            #                                       overwrite   = overwrite  )
        if _type == 'LSF':
            pass
            #self.build_peptide_from_sequence_AMBER(sequence    = sequence   , 
            #                                       force_field = force_field, 
            #                                       overwrite   = overwrite  )
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

    
    def import_CMAP (self, cmap = None, log = False, p_cmap = False):
        """ Function doc """
        
        self.energy_model_parameters['R_contact'] = cmap.cutoff
        self.R_contact                            = cmap.cutoff
        _type                                     = cmap._type 
        cmap                                      = cmap.cmap 
        
        #cmap.number_of_contacts = 0
        
        if len(cmap) != len(self.residues):
            if log:
                print 'error - cmap wrong size.' 
                print 'cmap size  : ', len(cmap) 
                print 'system size: ', len(self.residues) 
        else:
            if log:
                print '\n' 
                print 'cmap size  : ', len(cmap) 
                print 'system size: ', len(self.residues)
            
            self.cmap = cmap
            
            if p_cmap:
                for index in range(0, len(cmap)):
                    print cmap[index]


    def import_fixed_from_string(self, fixed=None):
        """ Function doc """

        #-------------------importing from porter------------------------------
        n = 0
        for aa in fixed:
            if aa != '0':
                self.fixed_residues.append(n)
            n += 1
        #----------------------------------------------------------------------

    def Status (self, parameters = True, coordinates = True ):
        """ Function doc """
        #print 'Bond_Stretching  = ', self.Bond_Stretching
        #print 'Angle_Bending    = ', self.Angle_Bending
        #print 'Improper_Torsion = ', self.Improper_Torsion
        #print 'Torsional_Angle  = ', self.Torsional_Angle
        #print 'Van_der_Waals    = ', self.Van_der_Waals
        #print 'Charge_Charge    = ', self.Charge_Charge
        #print 'Total_Energy     = ', self.Total_Energy

        #print 'Fixed_residues   = ', self.fixed_residues
        n_atoms = 0
        for res in self.residues:
            for atom in res.atoms:
                #print atom.name, atom.pos
                n_atoms += 1
        
        text = '''
--------------------------------------------------------------------------------
                        Summary for System "%s"
--------------------------------------------------------------------------------

---------------------------- Atom Container Summary ----------------------------
Number of Residues     =  %10d   Number of Atoms      =  %10d
--------------------------------------------------------------------------------

        '''%(self.name, len(self.residues), n_atoms)
        
        


        #for item in self.energy_model_parameters:
        #    print '%10s = %10.5f' %(item, self.energy_model_parameters[item])
        
        
        
        
        text += '''
-------------------------------------------------------------------------------       
                        Summary for Energy Model "%s"
-------------------------------------------------------------------------------       
ANGLE          =  %19.10f    ESURF           =  %19.10f
DIHED          =  %19.10f    NB              =  %19.10f
EEL            =  %19.10f    VDWAALS         =  %19.10f
EELEC          =  %19.10f    EGB'            =  %19.10f 
-------------------------------------------------------------------------------       
        ''' % ( self.energy_model, 
            self.energy_model_parameters['ANGLE'  ],
            self.energy_model_parameters['ESURF'  ],
            self.energy_model_parameters['DIHED'  ],

            self.energy_model_parameters['NB'     ],

            self.energy_model_parameters['EEL'    ],
            self.energy_model_parameters['VDWAALS'],

            self.energy_model_parameters['EELEC'  ],
            self.energy_model_parameters['EGB'    ],
            #self.energy_model_parameters['CONTACT'],
            #self.energy_model_parameters['AB'     ]
               ) 
      

        text += '''
---------------------- AMBER Single Point Calculations ------------------------
IGB                              = % 8d  Cutoff            = %8.2f
Saltcon                          = % 8.5f  GBSA              = %8d
rgbmax                           = % 8.2f  Surftern          = %8.4f
------------------------------------------------------------------------------- 
        '''% (
               self.energy_model_parameters['igb'    ],
               self.energy_model_parameters['cut'    ],
               self.energy_model_parameters['saltcon'],
               self.energy_model_parameters['gbsa'   ],
               self.energy_model_parameters['rgbmax' ],
               self.energy_model_parameters['surften'])

       
        
        text += '''
--------------------------- Calpha Calculations -------------------------------
AB weight                   = %8.5f  Cutoff            = %8.2f
Hydrofobic type             = %8s  
------------------------------------------------------------------------------- 
        '''% (
               self.energy_model_parameters['AB'      ],
               self.energy_model_parameters['R_cutoff'],
               'simple' , )
               
        
        text += '''
-------------------- Contact Energy Model Calculations ------------------------
CONTACT weight              = %8.5f  Cutoff            = %8.2f
Matrix type                 = %8s  
------------------------------------------------------------------------------- 
        '''% (
               self.energy_model_parameters['CONTACT'      ],
               self.energy_model_parameters['R_contact'],
               None , )


       
        print text 


    def import_SS_from_string(self, ss = None):
        # TIDQWLLKNAKEDAIAELKKAGITSDFYFNAINKAKTVEEVNALKNEILKAHA
        # CCHHHHHHHHHHHHHHHHHHCCCCCHHHHHHHHHCCCHHHHHHHHHHHHHHCC
        if len(self.residues) != len(ss):
            return False
        
        else:
           
            for i in range (0,len(self.residues)):
                if ss[i] == 'C':
                    pass
                
                if ss[i] == 'H':
                    from Geometry                 import *
                    phi_final_angle = set_phi_psi_dihedral( molecule=self, resi=i, bond='PHI',angle = -57 )
                    psi_final_angle = set_phi_psi_dihedral( molecule=self, resi=i, bond='PSI',angle = -47 )
                else:
                    pass

        
    def import_restraints_from_porter_file (self, filein= None):
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

