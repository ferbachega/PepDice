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



class Residue(object):
    """ Class doc """

    def __init__(self, id=0,  name='UNK'):
        """ Class initialiser """
        self.id = id
        self.atoms = []
        self.bonds = []
        self.name = name
        
        self.CHI1 = None
        self.CHI2 = None
        self.CHI3 = None
        self.CHI4 = None
        self.CHI5 = None
        
        # PHI
        self.phi_restraint_angle  = 0.0
        self.phi_restraint_weight = 0.0
        
        # PSI
        self.psi_restraint_angle  = 0.0
        self.psi_restraint_weight = 0.0

        # OMEGA
        self.omega_restraint_angle  = 0.0
        self.omega_restraint_weight = 0.0

        # CHI's
        self.chi1_restraint_angle  = 0.0
        self.chi2_restraint_angle  = 0.0
        self.chi3_restraint_angle  = 0.0
        self.chi4_restraint_angle  = 0.0
        self.chi5_restraint_angle  = 0.0
        self.chi1_restraint_weight = 0.0
        self.chi2_restraint_weight = 0.0
        self.chi3_restraint_weight = 0.0
        self.chi4_restraint_weight = 0.0
        self.chi5_restraint_weight = 0.0

        #self.ss_restraint   = [None, None, None]
        #self.w_ss_restraint = [0.0 , 0.0,  0.0 ]
     
     
     
     
     
     
     
        
'''
secondary_structure

    enabled = False Turn on secondary structure restraints (main switch)
    protein
        enabled = True Turn on secondary structure restraints for protein
        search_method = *ksdssp mmtbx_dssp from_ca cablam Particular method to search protein secondary structure.
        distance_ideal_n_o = 2.9 Target length for N-O hydrogen bond
        distance_cut_n_o = 3.5 Hydrogen bond with length exceeding this value will not be established
        remove_outliers = True If true, h-bonds exceeding distance_cut_n_o length will not be established
        helix
            serial_number = None
            helix_identifier = None
            enabled = True Restrain this particular helix
            selection = None
            helix_type = *alpha pi 3_10 unknown Type of helix, defaults to alpha. Only alpha, pi, and 3_10 helices are used for hydrogen-bond restraints.
            sigma = 0.05
            slack = 0
            top_out = False
            hbond
                donor = None
                acceptor = None 
        sheet
            enabled = True Restrain this particular sheet
            first_strand = None
            sheet_id = None
            sigma = 0.05
            slack = 0
            top_out = False
            strand
                selection = None
                sense = parallel antiparallel *unknown
                bond_start_current = None
                bond_start_previous = None 
            hbond
                donor = None
                acceptor = None 
'''
