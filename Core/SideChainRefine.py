#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  SideChainRefine.py
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

#----------------------------------------------------#
from Geometry    import *                            #
#----------------------------------------------------#
#----------------------------------------------------#
from multiprocessing import Pool                     #
#----------------------------------------------------#



def compute_energy(parameters):
    system = parameters['system']
    resi   = parameters['resi']
    rotamer= parameters['rotamer']
    pn     = parameters['pn']
    
    results= {
             }

    set_side_chain_rotamer(molecule=system, resi=resi, rotamer=rotamer)
    energy            = system.energy(pn = pn)
    
    #save_XYZ_to_file(system, TRAJECTORY)
    results[energy] = system.get_coordinates_from_system()
    #results['coords'] = system.get_coordinates_from_system()
    return results


def optimize_side_chain(system = None, resi = None, rotamer_list = None, pool = 8):
    """ Function doc """
    if rotamer_list == None:
        from AATorsions  import ROTAMER_LIST 
    else:
        ROTAMER_LIST = rotamer_list
        
    name    = system.residues[resi].name
    res     = system.residues[resi]
    process = []
    
    print name
    
    if name in ['HSD', 'HSE', 'HDP', 'HIE', 'HID']:
        name = 'HIS'
    #res_rotamers =  ROTAMER_LIST[name]
    pn = 1
    for key in ROTAMER_LIST[name]:
        parameters = {}
        parameters['system']  = system
        parameters['resi']    = resi
        parameters['rotamer'] = ROTAMER_LIST[name][key]
        parameters['pn']      = pn 
        process.append(parameters)
        pn += 1

    #print len(process)                 
    p = Pool(pool)                      
    
    results = p.map(compute_energy, process)     
    p.close()
    p.join()

    RESULTS = {}
    for result in results:          
        a = result.keys()
        if a[0] != None:
            if result[a[0]] == None:
                pass
            else:
                RESULTS[a[0]] = result[a[0]]
    return RESULTS

