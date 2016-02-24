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
        
    with open(filename, "r") as crdin:
        text  = crdin.read()
        text.replace('\n', '')
        text  = text.split()
        
        pos = []
        n = 0
        j = 1
        
        for coord in text[2:(len(text)+1)]:
           pos.append(coord)
           n = n + 1
           #print coord
           
           if n == 3:
                j   = j + 1
                newcoords.append(pos)
                n   = 0
                pos = []
        
        #print 'jota', j , len(newcoords), newcoords[-1], text[-1]
    
    n = 0
    newcoords
    for residue in molecule.residues:
        for atom in residue.atoms:
            #print n, atom.id, newcoords[n][0], newcoords[n][1], newcoords[n][2]
            
            atom.pos[0] = newcoords[n][0]
            atom.pos[1] = newcoords[n][1]
            atom.pos[2] = newcoords[n][2]
            n = n+1

        

        
        
        #for line in crdin:
        #    line2 = line.split()
        #    if len(line2) >= 3:
        #        print 
        

        #for residue_i in molecule.residues:
        #    for atom_i in residue_i.atoms:
        #        n = n +1
        #        x = float(atom_i.pos[0])
        #        y = float(atom_i.pos[1])
        #        z = float(atom_i.pos[2])
        #        
        #        if n == 2:
        #        
        #            text += "%12.7f%12.7f%12.7f\n" %(x ,y,z)
        #        
        #            n = 0
        #        
        #        else:
        #            text += "%12.7f%12.7f%12.7f" %(x ,y,z)
        #
        #output_file.write(text)
        #output_file.close()
        #pass

#a  = 'b'
#load_CRD_from_file (filename = '/home/fernando/programs/PepDice/Examples/1GAB/amber/1gab_amber.crd')


