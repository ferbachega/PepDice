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




class Atom(object):
    """class to store info about an atom"""

    def __init__(self, id=0,
                 a_number=None,
                 name=None,
                 sigma=0.3405,
                 epsilon=0.1,  # 119.8 ,
                 charge=1.0,
                 pos=[0, 0, 0]):

        # 119.8 K and  0.3405 nm)

        # PDB data
        self.id = id  # index
        self.name = name  # nome do atomo eg  "OH" ou "CA"
        # self.fname          = None          #  nome completo do atomo " oxygen"
        # self.resn           = resn          #  nome do residuo
        # self.resi           = resi          #  numero do residuo
        self.element = None

        # FF paramters
        self.charge = charge  # carga
        self.sigma = sigma
        self.epsilon = epsilon
        self.ff_atom_number = None
        self.ff_atom_name = None
        self.atomic_number = None
        self.mass = None

        # coordinates and energy
        self.pos = pos
        self.actual_energy = 0
