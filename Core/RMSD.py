#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  RMSD.py
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


def compute_RMSD (crd1 = None, crd2= None):
    """ Function doc """
    N    = len(crd1)
    soma = 0.0
    for i in range(0, N):
        crd1_x = crd1[i][0]
        crd1_y = crd1[i][1]
        crd1_z = crd1[i][2]
        
        crd2_x = crd2[i][0]
        crd2_y = crd2[i][1]
        crd2_z = crd2[i][2]
        #print crd1[i], crd2[i] 
        soma += (crd1_x - crd2_x)**2 + (crd1_y - crd2_y)**2 + (crd1_z - crd2_z)**2
        #print soma
     
    #print soma/N
    rmsd = (soma/N)**0.5
    #print 'rmsd = ', rmsd
    return rmsd
