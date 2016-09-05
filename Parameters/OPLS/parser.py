#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  parser.py
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


filein = open('oplsua.prm', 'r')


atom    = {}
vdw     = {}
bond    = {}
angle   = {}
torsion = {}
imptors = {}

for line in filein:
    line2 = line.split()
  
    if len(line2) > 0 :
        if line2[0] == 'atom':
            atom[line2[1]] = [line2[2], line2[-3], line2[-2], line2[-1]]
            
        if line2[0] == 'vdw':
            vdw[line2[1]] = [line2[2],line2[3]]
            
        if line2[0] == 'bond':
            bond[(line2[1],line2[2])] = [line2[3],line2[4]]

        
        if line2[0] == 'angle':
            angle[ (line2[1],line2[2],line2[3]) ] = [line2[4],line2[5]]

        
        
        if line2[0] == 'torsion':
            try:
                torsion[ (line2[1],line2[2],line2[3],line2[4]) ] = [line2[5],line2[6]]
            except:
                print (line2[1],line2[2],line2[3],line2[4]), 'has no parameters'

        if line2[0] == 'imptors':
            try:
                imptors[ (line2[1],line2[2],line2[3],line2[4]) ] = [line2[5],line2[6],line2[7]]
            except:
                print (line2[1],line2[2],line2[3],line2[4]), 'has no parameters'

#for i in atom:
#    print i, atom[i], vdw[i]
#    
#for i in bond:
#    print i, bond[i]
#
for i in angle:
    print i, atom[i[0]][0], atom[i[1]][0],atom[i[2]][0]
