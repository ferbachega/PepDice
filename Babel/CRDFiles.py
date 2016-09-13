#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  CRDFiles.py
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


def save_CRD_to_file(molecule, filename):
    with open(filename, "w") as output_file:
        text = '\n'
        text += ' ' + str( molecule.residues[-1].atoms[-1].id) + '\n'
        n = 0
        for residue_i in molecule.residues:
            for atom_i in residue_i.atoms:
                n = n +1
                x = float(atom_i.pos[0])
                y = float(atom_i.pos[1])
                z = float(atom_i.pos[2])
                
                if n == 2:
                
                    text += "%12.7f%12.7f%12.7f\n" %(x ,y,z)
                
                    n = 0
                
                else:
                    text += "%12.7f%12.7f%12.7f" %(x ,y,z)

        output_file.write(text)
        output_file.close()
        pass



def load_CRD_from_file (molecule = None, filename = None):
    newcoords  =  []
    
    '''   3.3380001   2.4640000  -0.2370000   3.7279999   1.7250000  -0.0340000'''    
    
    with open(filename, "r") as crdin:
        text  = crdin.readlines()
        #text  = text.readlines() 
        pos = []
        n   = 0
        j   = 1
        
        for line in text[2:]:
            line.replace('\n', '')
            
            #print line[0:12], line[12:24], line[24:36], line[36:48], line[48:60], line[60:72]

            pos.append(line[0:12] )
            pos.append(line[12:24])
            pos.append(line[24:36])
            pos.append(line[36:48])
            pos.append(line[48:60])       
            pos.append(line[60:72])
            
    n = 0
    if molecule != None:
        for residue in molecule.residues:
            for atom in residue.atoms:
                
                atom.pos[0] = pos[n]
                n += 1
                
                atom.pos[1] = pos[n]
                n += 1
                
                atom.pos[2] = pos[n]
                n += 1




load_CRD_from_file(molecule = None, filename = '/home/fernando/programs/pepdice/Examples/SinglePoint1.crd')
