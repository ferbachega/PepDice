#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  MonteCarlo.py
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
#-------------------------------------------
import math
#-------------------------------------------




def metropolis_acceptance_test (energy = None        ,
                       previous_energy = None        ,
                           temperature = 273.15      ,
                                    Kb = 1           ,
                                random = None): #0.0019872041):
    """ Function doc """
    if energy < previous_energy:
        return True
    
    else:
        dE = (energy - previous_energy)
        Px = math.exp(-1 * dE / (temperature))
        X  = random.uniform(0, 1)
        return X <= Px



    




