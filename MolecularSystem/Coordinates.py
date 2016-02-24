#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  Coordinates.py
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

class Coordinates:
    """ Class doc """
    
    def __init__ (self):
        """ Class initialiser """
        pass
        
    def import_coordinates_to_system(self, coordinates):
        """ Function doc """
        # for key in self.atoms.keys():
        for residue_i in self.residues:
            for atom_i in residue_i.atoms:
                index_i = atom_i.id
                atom_i.pos[0] = coordinates[index_i - 1][0]
                atom_i.pos[1] = coordinates[index_i - 1][1]
                atom_i.pos[2] = coordinates[index_i - 1][2]

    def get_coordinates_from_system(self):  # molecule = None):
        """ Function doc """
        # pprint(ff.biotypes)
        initial_coordinates = []
        for residue_i in self.residues:
            for atom_i in residue_i.atoms:
                initial_coordinates.append(
                    [atom_i.pos[0], atom_i.pos[1], atom_i.pos[2]])

        return initial_coordinates
