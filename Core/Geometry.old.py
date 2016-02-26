#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  Geometry.py
#
#  Copyright 2015 Fernando Bachega <fernando@Fenrir>
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

#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  geometry.py
#
#  Copyright 2015 farminf <farminf@farminf-3>
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

import numpy as np
import math


from Bio.PDB import PDBParser
#from force_field import ForceField
from pprint import pprint

import random
import numpy as np
import math
from Vectors import *
#from MolecularSystem import *
#from PDBFileReader   import *
#from XYZFileWriter   import *



def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation 
    about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    theta = np.asarray(theta)
    axis = axis / math.sqrt(np.dot(axis, axis))
    a = math.cos(theta / 2)
    b, c, d = -axis * math.sin(theta / 2)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                     [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                     [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])

def center_atom(molecule, x, y, z):

    for residue_i in molecule.residues:
        for atom_i in residue_i.atoms:
            atom_i.pos[0] = atom_i.pos[0] - x
            atom_i.pos[1] = atom_i.pos[1] - y
            atom_i.pos[2] = atom_i.pos[2] - z

    return molecule

def move_atom(molecule=None, index=None, delta_x=0, delta_y=0, delta_z=0):
    """ Function doc """
    atom = molecule.atoms[index]
    atom.pos[0] = atom.pos[0] + delta_x
    atom.pos[1] = atom.pos[1] + delta_y
    atom.pos[2] = atom.pos[2] + delta_z

def Compute_Rab(atom_a, atom_b):
    """ Function doc 
    Calcula o vetor R que eh a subtracao: vetor a  - vetor b
    """
    x_a = atom_a.pos[0]
    x_b = atom_b.pos[0]

    y_a = atom_a.pos[1]
    y_b = atom_b.pos[1]

    z_a = atom_a.pos[2]
    z_b = atom_b.pos[2]

    x = x_a - x_b
    # if x >= cutoff:
    #    return False
    y = y_a - y_b
    # if y >= cutoff:
    #    return False
    z = z_a - z_b
    # if z >= cutoff:
    #    return False
    return [x, y, z]





'''
----------------------------------------------------
                   SIDE CHAIN
----------------------------------------------------
'''
def computeCHI (molecule=None, resi=1, bond='CHI1'):
    """ Function doc """
   
    sideChain = []
    residue   = molecule.residues[resi]
    torsions  = molecule.torsions[residue.name]

    try:
        CHI = molecule.torsions[residue.name][bond]
        
        for atom in residue.atoms:
            if atom.name == CHI[0]:
                #print atom.name
                a1 = atom
            if atom.name == CHI[1]:
                #print atom.name
                a2 = atom
            if atom.name == CHI[2]:
                #print atom.name
                a3 = atom
            if atom.name == CHI[3]:
                #print atom.name
                a4 = atom 
        
        angle  =   dihedral(
                            a1.pos ,
                            a2.pos ,
                            a3.pos ,
                            a4.pos)
        
        return angle *57.324840764
    except:
        pass

def rotate_side_chain (molecule=None, resi=1, bond='CHI1', theta=0.017444445, steps=1):
    """ 
    This function rotates the side chain of a given amino acid , theta = 0.017444445 means 1o (um grau)
    to rotate 10o -> theta*10
    """
    index =  resi

    if resi in molecule.fixed_residues:
        return 0
        
    sideChain = []
     
    res      = molecule.residues[resi]
    torsions = molecule.torsions[res.name]
    
    if bond in molecule.torsions[res.name]:
        CHI      = molecule.torsions[res.name][bond]

        for atom in res.atoms:
            if atom.name == CHI[0]:
                #print atom.name
                a1 = atom
            if atom.name == CHI[1]:
                #print atom.name
                a2 = atom
            if atom.name == CHI[2]:
                #print atom.name
                a3 = atom
            if atom.name == CHI[3]:
                #print atom.name
                a4 = atom 

        angle  =   dihedral(
                            a1.pos ,
                            a2.pos ,
                            a3.pos ,
                            a4.pos)
        
        subcoord = a2.pos
        x = subcoord[0]
        y = subcoord[1]
        z = subcoord[2]
        molecule = center_atom(molecule=molecule, x=x, y=y, z=z)
        axis = a3.pos
        
        for atom_i in res.atoms:
            if atom_i.name in molecule.FIX_atoms_CHI[bond]:
                pass
            else:
                atom_i.pos = np.dot(rotation_matrix(axis,theta),
                                        atom_i.pos
                                        )   

def set_chi_dihedral (molecule=None, resi=1, bond='CHI1', angle = 0.0 ):
    """ 
    This function changes a CHI angle to a specific value provided by the user.
    """
    initial_angle = computeCHI(molecule = molecule, 
                               resi     = resi, 
                               bond     = bond)
    
    #print initial_angle
    delta_angle = (angle - initial_angle) #* -1
    theta       = delta_angle* 0.017444445 #  theta equivalente a 1 grau
    
    
    rotate_side_chain (molecule = molecule  , 
                       resi     = resi      , 
                       bond     = bond      , 
                       theta    = theta     , 
                       steps    = 1)
    
    final_angle =   computeCHI(molecule = molecule, 
                               resi     = resi, 
                               bond     = bond)
    
    print 'initial :', initial_angle, 'final', final_angle
    return final_angle
    
def set_side_chain_rotamer (molecule=None, resi=1, rotamer=None):
    '''
    This function assigns a specific rotamer for a given amino acid.
    '''
    for chi in rotamer:
        print chi, rotamer[chi]
        
        if rotamer[chi] == None:
            pass
        else:
            #print rotamer
            if chi == 'THETA':
                pass
            else:
                angle = rotamer[chi]
                bond  =  chi
                set_chi_dihedral (molecule = molecule    , 
                                      resi = resi        , 
                                      bond = chi         , 
                                     angle = rotamer[chi])
        





   
   

def rotate_Calpha_dihedral(molecule, axis, theta, window):
    """ Function doc """
    for residue_i in molecule.residues:
        for atom_i in residue_i.atoms:

            if atom_i.id in window:
                atom_i.pos = np.dot(rotation_matrix(axis,
                                                    theta),
                                    atom_i.pos
                                    )
    
def rotate_backbone(molecule=None, resi=1, bond='PSI', theta=0, steps=1):
    """ Function doc """
    
    index =  resi
    if resi in molecule.fixed_residues:
        return 0

   
    sideChain = []
    residue = molecule.residues[resi]
        
    for atom_i in residue.atoms:

        if atom_i.name == 'CA':
            CA = atom_i

        elif atom_i.name == 'N':
            N = atom_i

        elif atom_i.name in ['H', 'HT1', 'HT2', 'HT3', 'HN']:
            H = atom_i

        elif atom_i.name == 'C':
            C = atom_i

        elif atom_i.name == 'O':
            O = atom_i
            
        else:
            sideChain.append(atom_i.id)
            
    if bond == 'PSI':
        subcoord = CA.pos
        x = subcoord[0]
        y = subcoord[1]
        z = subcoord[2]

        molecule = center_atom(molecule=molecule, x=x, y=y, z=z)
        axis = C.pos
        window = range(0, CA.id + 1)

        for i in sideChain:
            window.append(i)

        for i in range(0, steps):
            rotate_Calpha_dihedral(molecule, axis, theta, window)

    if bond == 'PHI':
        subcoord = CA.pos  # N - nitrogenio
        x = subcoord[0]
        y = subcoord[1]
        z = subcoord[2]
        molecule = center_atom(molecule=molecule, x=x, y=y, z=z)

        axis = N.pos     # CA - C alpha
        window = range(0, N.id + 1)
        try:
            window.append(H.id)
        except:
            pass

        for i in range(0, steps):
            rotate_Calpha_dihedral(molecule, axis, theta, window=window)
    
def set_phi_psi_dihedral (molecule=None, resi=1, bond='PSI', angle = 0.0 ):
    """ Function doc """
    
    #print resi, molecule.residues[resi].name 
    #print bond
    #print angle
    
    
    initial_angle =  computePhiPsi(molecule = molecule, 
                                   resi     = resi, 
                                   bond      = bond)
    if initial_angle == 0.0:
        initial_angle = 180.0

    if initial_angle == None:
        #print bond, 'False'
        return False

    if bond == 'PSI':
        delta_angle = (angle -  initial_angle) * -1
        theta       = delta_angle* 0.017444445 #  theta equivalente a 1 grau 
    else:
        delta_angle = (angle -  initial_angle) 
        theta       = delta_angle* 0.017444445 #  theta equivalente a 1 grau 
        
    rotate_backbone(molecule = molecule   , 
                    resi     = resi      , 
                    bond     = bond       , 
                    theta    = theta, 
                    steps    = 1)
    
    final_angle = computePhiPsi     (molecule  = molecule, 
                                      resi     = resi, 
                                      bond     = bond)
    
    #print bond, resi, molecule.residues[resi].name, initial_angle, final_angle, angle, delta_angle
    return final_angle
    











def computePhiPsi (molecule=None, resi=1, bond='PSI'):
    """ Function doc """
    C1  = None
    N2  = None
    CA2 = None
    C2  = None
    N3  = None
   
    
    # obtaining the CA N C positions residue n 
    residue = molecule.residues[resi]
    for atom in residue.atoms:
        if atom.name == 'CA':
            CA2 = atom

        if atom.name == 'N':
            N2 = atom
            
        if atom.name == 'C':
            C2 = atom
    #print 'CA', CA2, 'N2', N2, 'C2', C2
    
    #----------------------------------------------
    
    # obtaining the C positions residue n - 1 
    if resi == 0:                             
        phi = False
        C1  = None
        CA1 = None
        pass
            
    else:                                     
        residue = molecule.residues[resi - 1] 
        for atom in residue.atoms:            
            if atom.name == 'C':              
                C1 = atom                     
                                              
            if atom.name == 'CA':             
                CA1 = atom                    
        phi = True
    #----------------------------------------------
    
    #print 'CA1', CA1,'C1', C1, 'psi', phi



    #----------------------------------------------
    # obtaining the C positions residue n + 1 
    try:
        residue = molecule.residues[resi + 1]
        for atom in residue.atoms:
            if atom.name == 'N':
                N3 = atom
            if atom.name == 'CA':
                CA3 = atom
        
        psi = True
        ome = True
    except:
        psi = False
        ome = False
    #----------------------------------------------



    
    #print phi, psi, bond
    if phi:
        if bond == 'PHI':
            angle = dihedral(C1.pos  ,
                             N2.pos  ,
                             CA2.pos ,
                             C2.pos  )
            #print 'PHI', C1.name,N2.name,CA2.name,C2.name, angle*57.324840764
            return angle*57.324840764
    if psi:
        if bond == 'PSI':
            angle = dihedral(
                             N2.pos  ,
                             CA2.pos ,
                             C2.pos  ,
                             N3.pos)
            #print 'PSI', N2.name,CA2.name,C2.name,N3.name, angle*57.324840764
            return angle*57.324840764
        
    if ome:
        if bond == 'OMEGA':
            angle = dihedral(
                             CA1.pos  ,
                             C1.pos ,
                             N2.pos  ,
                             CA2.pos)
            return angle*57.324840764




























'''
if __name__ == '__main__':
    """Code to test the atom, bond and molecule classes"""

    m = load_CHARMM_molecular_system(pdb='examples/HcH_autopsf.pdb',
                                     psf='examples/HcH_autopsf.psf',
                                     param='ff/charmm/par_all27_prot_lipid_na.inp')

    save_XYZ_to_file(m, 'Geometry_test.xyz')
    for i in range(0, 100):
        rotate_backbone(molecule=m, residue=12,
                        bond='psi', theta=0.01, steps=1)
        save_XYZ_to_file(m, 'Geometry_test.xyz')
    for i in range(0, 100):
        rotate_backbone(molecule=m, residue=12,
                        bond='phi', theta=0.01, steps=1)
        save_XYZ_to_file(m, 'Geometry_test.xyz')
    #rotate_bond (molecule = m, residue = 2, bond = 'psi', theta = 0.01 , steps = 100)
    # print initial_coordinates



def center_atom (molecule, x, y, z):

	
	for residue in molecule.Residues:
		#print residue.Name
		for atom in residue.Atoms:
	
			#print atom.Name, atom.coordinates#, subcoord
			atom.coordinates[0] = atom.coordinates[0] - x
			atom.coordinates[1] = atom.coordinates[1] - y
			atom.coordinates[2] = atom.coordinates[2] - z
			#print atom.Name, atom.coordinates#, subcoord
	
	
	return molecule
	
def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation 
    about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    theta = np.asarray(theta)
    axis = axis/math.sqrt(np.dot(axis, axis))
    a = math.cos(theta/2)
    b, c, d = -axis*math.sin(theta/2)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])




def move_fragment (molecule, delta, window):
	""" Function doc """
	
	for residue in molecule.Residues[window[0]:window[1]]:
		for atom in residue.Atoms:
			atom.coordinates += delta
	
	

def rotate_Calpha_dihedral (molecule, axis, theta, window):
    """ Function doc """
    #v = [3, 5, 0]
    #axis = [4, 4, 1]
    #theta = 1.2 

    #print(np.dot(rotation_matrix(axis,theta), v)) 
    for residue in molecule.Residues:
        for atom in residue.Atoms:
            if atom.id in window:
                atom.coordinates = np.dot(
                    rotation_matrix(axis,theta), 
                    atom.coordinates
                )
			
			
def main():
	
	return 0

if __name__ == '__main__':
	main()
'''