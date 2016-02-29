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
from pprint import pprint
import time

import random
import numpy as np
import math
import os
from Geometry import *
from XYZFiles import save_XYZ_to_file






def MC_test_energy (energy = None        , 
           previous_energy = None        ,
               temperature = 273.15      ,
                        Kb = 0.0083144621):
    """ Function doc """
    if energy < previous_energy:
        return True
    else:
        DG = (energy - previous_energy)
        Px = math.exp(-1 * DG / (Kb * temperature))
        X  = random.uniform(0, 1)
        return X <= Px
            

def rotate_backbone_attempt (molecule = None, 
                                 resi = None, 
                                 bond = None,
                                theta = None,  
                      previous_energy = None,
                          temperature = None,
                                   pn = None):
    
    rotate_backbone(molecule = molecule, 
                        resi = resi    ,
                        bond = bond    , 
                       theta = theta   )

    energy = molecule.energy(pn =pn)
    
    if energy:
        if MC_test_energy (energy = energy         , 
                 previous_energy = previous_energy ,
                     temperature = 273.15          ,
                              Kb = 0.0083144621    ):    
            return energy

        else:
            return False
    else:
        return False


def insert_fragment (molecule = None, fragment = None, sidechain = False):
    """ Function doc """
    for key in fragment:
        #PSI = fragment[key]['PSI']
        #PHI = fragment[key]['PHI']
        for bond in ['PSI','PHI']:
            #print  key, bond,  fragment[key][bond]
            set_phi_psi_dihedral (molecule = molecule, 
                                      resi = key     , 
                                      bond = bond    , 
                                      angle = fragment[key][bond])
        
            
        if sidechain:
            try:
                for bond in ['CHI1','CHI2','CHI3','CHI4','CHI5']:
                    if fragment[key][bond]:
                        set_chi_dihedral (molecule  = molecule, 
                                              resi  = key, 
                                              bond  = bond, 
                                              angle = fragment[key][bond])
            except:
                print 'failed sidechain'

def monte_carlo_dic (parameters):
    """ Function doc """
    results = monte_carlo(
                          molecule           = parameters['molecule'          ],
                          temperature        = parameters['temperature'       ],
                          Kb                 = parameters['Kb'                ],
                          angle_range        = parameters['angle_range'       ],
                          nSteps             = parameters['nSteps'            ],
                          fragment_rate      = parameters['fragment_rate'     ],
                          fragment_sidechain = parameters['fragment_sidechain'],
                          PhiPsi_rate        = parameters['PhiPsi_rate'       ],
                          trajectory         = parameters['trajectory'        ],
                          pn                 = parameters['pn'                ],
                          )
    return results


def monte_carlo(molecule           = None                       ,
                temperature        = 1000                       ,
                Kb                 = 0.0083144621               ,
                angle_range        = 1                          ,
                nSteps             = 10000                      ,
                fragment_rate      = 1.0                        , #between 0  and 1
                fragment_sidechain = False                      ,
                PhiPsi_rate        = 1.0                        ,
                trajectory         = 'MonteCarlo_trajectory.xyz',
                pn                 = 1                          ):


    n = 0
    previous_energy      = molecule.energy(pn = pn)
    previous_coordinates = molecule.get_coordinates_from_system()
    previous_fragment    = None
    
    for i in range(0, nSteps):
        #--------------------------------#
        #       attempted_accepted       #
        #--------------------------------#
        attempted_fragment = 0.0         #
        accepted_fragment  = 0.0         #
                                         #
        attempted_phi = 0.0              #
        accepted_phi  = 0.0              #
                                         #
        attempted_psi = 0.0              #
        accepted_psi  = 0.0              #
        #--------------------------------#
        
        
        fragment_acceptance = random.uniform(0, 1)
        #----------------------------------------------#
        #                  FRAGMENTS                   #
        #----------------------------------------------#
        # se o numero for menor ou igual a chance, entao um novo fragmento eh atribuido a estrutura
        if fragment_acceptance <= fragment_rate:
            attempted_fragment += 1
            fragment_index = random.randint(0, len(molecule.fragments)-1)            
            fragment = molecule.fragments[fragment_index]
            
            if fragment != previous_fragment:
                previous_fragment = fragment
                
                insert_fragment (molecule   = molecule, 
                                 fragment   = fragment,
                                 sidechain  = fragment_sidechain)
                                 
                energy = molecule.energy(pn = pn)
                if energy:
                    if MC_test_energy (energy = energy         , 
                               previous_energy = previous_energy ,
                                  temperature = temperature    ):
                                           
                        save_XYZ_to_file (molecule, trajectory)
                        previous_energy      = energy
                        previous_coordinates = molecule.get_coordinates_from_system()
                        accepted_fragment += 1
                        print "pn: {:<3d} step: {:5d} energy: {:<20.7f}fragment: {:<3d}".format(pn, i, energy , fragment_index)
                    else:
                        molecule.import_coordinates_to_system (previous_coordinates)
                        #print 'fragment: ',fragment_index, energy, len(fragment), 'failed'
                else:
                    molecule.import_coordinates_to_system (previous_coordinates)
        #print 'temp: = ', temperature, 'energy = ', previous_energy, 'acceptance ratio (phi) =', (accepted_fragment / attempted_fragment)

       
        #----------------------------------------------#
        #               PHI/PSI sampling               #
        #----------------------------------------------#
        PhiPsi_acceptance = random.uniform(0, 1)
        if PhiPsi_acceptance <= PhiPsi_rate:
            resi = random.randint(0, len(molecule.residues)-1)
            
            if resi in molecule.fixed_residues:
                pass
            else:
                #print resi
                for bond in ['PSI','PHI']:
                    
                    if bond == 'PSI':
                        attempted_phi += 1
                    if bond == 'PHI':
                        attempted_psi += 1
                    
                    #attempted_psi += 1
                    theta  = random.uniform(-1 * angle_range, angle_range)
                    theta = theta * 0.017444445
                    
                    #------------------------------------------#
                    #         rotate_backbone_attempt          #
                    #------------------------------------------#
                    
                    energy = rotate_backbone_attempt (molecule = molecule      , 
                                                        resi = resi            , 
                                                        bond = bond            ,
                                                       theta = theta           ,  
                                              previous_energy = previous_energy,
                                                 temperature = temperature     ,
                                                          pn = pn)
        
                    if energy:
                        save_XYZ_to_file (molecule, trajectory)
                        previous_energy      = energy
                        previous_coordinates = molecule.get_coordinates_from_system()
                        
                        if bond == 'PSI':
                            accepted_psi += 1
                        if bond == 'PHI':
                            accepted_phi += 1
                        print "pn: {:<3d} step: {:5d} energy: {:<20.7f}rotate_backbone theta: {:<6.3f}".format(pn, i, energy , theta*57.324)
                        #print '%5i %4i %10.4f %10.4f %3i' % (i, resi, theta*57.324, energy, temperature) #, previous_energy )  
                    else:
                        molecule.import_coordinates_to_system (previous_coordinates)
    return {'pn':pn, 'energy': previous_energy, 'coords': previous_coordinates, 'temperature': temperature }


def MC_replica_exchange (replicas = [], cpus = 8, Kb = 0.0083144621):
    from multiprocessing import Pool
    p = Pool(cpus)
    results = p.map(monte_carlo_dic, replicas)
    
    
    #--------------------------- Exchange ---------------------------#
    REPLICAS = {} 
    for result in results:
        print 'replica: %3i energy: %10.7f' %(result['pn'], result['energy'])#, len(result[2])
        REPLICAS[result['pn']] = {
                               'energy'     : result['energy'],
                               'coords'     : result['coords'],
                               'temperature': result['temperature']
                               }
    #----------------------------------------------------------------#

    for i in REPLICAS:
        for j in REPLICAS:
            
            if j == i:
                pass
            else:
                deltaG = REPLICAS[j]['energy'] - REPLICAS[i]['energy']
                div    = ((1/Kb *REPLICAS[j]['temperature']) - (1/Kb *REPLICAS[i]['temperature'])) 
                #print i, j , 'div', div * deltaG, (REPLICAS[i]['energy'] - REPLICAS[j]['energy']) *((1/Kb *REPLICAS[i]['temperature']) - (1/Kb *REPLICAS[j]['temperature'])) 
                #p      = math.exp(deltaG * div)
                #print i , j, deltaG, p
                
                
    
    #pprint(REPLICAS)

def run_MC_replica_exchange (
                            molecule           = None         ,
                            N_replicas         = 1            , # >= number of CPUs
                            min_temp           = 50           ,
                            max_temp           = 250          ,
                            PhiPsi_rate        = 0.1          , 
                            max_angle_range    = 5            ,
                            trajectory         = 'MC_replica_',      
                            Kb                 = 0.0083144621 ,
                            nSteps             = 1000         ,
                            fragment_rate      = 0.5          ,
                            fragment_sidechain = True         ,
                            ):
    """ Function doc """
    
    temperature_factor = (max_temp-min_temp)/N_replicas
    
    
    replicas   = []
    for i in range(1, N_replicas + 1):
        try:
            os.remove(TRAJECTORY +str(i)+'.xyz' )
        except:
            pass
        parameters = {} 
        parameters['molecule'          ] = molecule      
        parameters['temperature'       ] = min_temp       
        parameters['Kb'                ] = Kb
        parameters['nSteps'            ] = nSteps     
        parameters['fragment_rate'     ] = fragment_rate        
        parameters['fragment_sidechain'] = fragment_sidechain        
        parameters['PhiPsi_rate'       ] = PhiPsi_rate       
        parameters['angle_range'       ] = max_angle_range          
        parameters['trajectory'        ] = trajectory +str(i)+'.xyz' 
        parameters['pn'                ] = i
        replicas.append(parameters)
        
        min_temp += temperature_factor
        
    print 'number of residues: ' , len(molecule.residues)
    print 'number of fragments:' , len(molecule.fragments)
    print 'number of replicas: ' , len(replicas)
    #monte_carlo_dic(replicas[0])
    MC_replica_exchange(replicas= replicas, cpus = N_replicas)









def monte_carlo_side_chain (molecule  = None,
                           initial_T  = 1000,
                           final_T    = 1   ,
                           angle_range= 20  ,
                           nSteps     = 1000,
                           trajectory='MonteCarlo_trajectory.xyz'):
    
    Kb              = 0.0083144621
    temp            = initial_T
    gamma           = -1 * (math.log(float(final_T) / float(initial_T))) / nSteps
    tau_angle_range = (float(angle_range) / nSteps) * -1
    n = 0
    
    initial_energy = molecule.energy()
    initial_coordinates = molecule.get_coordinates_from_system()
    print 'gamma = ', gamma

    for i in range(0, nSteps):
        
        for i in range(0, len(molecule.residues)):
            name         = system.residues[i].name
            res          = system.residues[i]
            if name in ['HSD', 'HSE', 'HDP', 'HIE', 'HID']:
                name = 'HIS'
            res_rotamers =  ROTAMER_LIST[name]
            
            for key in res_rotamers:
                set_side_chain_rotamer(molecule=system, resi=i, rotamer=res_rotamers[key])
                
                energy = system.energy()
                if energy:   
                    if energy <= initial_energy:
                        save_XYZ_to_file(molecule, trajectory)
                        initial_energy = energy
                        initial_coordinates = molecule.get_coordinates_from_system()
                    else:
                        DG = (energy - initial_energy)
                        Px = math.exp(-1 * DG / (Kb * temp))
                        X = random.uniform(0, 1)
                        if X <= Px:
                            save_XYZ_to_file(molecule, trajectory)
                            initial_energy = energy
                            initial_coordinates = molecule.get_coordinates_from_system()
                        else:
                            molecule.import_coordinates_to_system(
                                initial_coordinates)
                else:
                    print 'Clash!'
        print 'temp: = ', temp, 'energy = ', initial_energy
        temp = float(initial_T) * math.exp(-1 * gamma * i)

'''
def monte_carlo(molecule   = None,
               initial_T   = 1000,
               final_T     = 1   ,
               angle_range = 20  ,
               nSteps      = 1000,
               trajectory='MonteCarlo_trajectory.xyz'):
    """ Function doc """
    #
    Kb = 0.0083144621
    temp = initial_T

    gamma = -1 * (math.log(float(final_T) / float(initial_T))) / nSteps
    tau_angle_range = (float(angle_range) / nSteps) * -1

    n = 0
    initial_energy = molecule.energy()
    initial_coordinates = molecule.get_coordinates_from_system()
    print 'gamma = ', gamma

    
    for i in range(0, nSteps):
        attempted_phi = 0.0
        accepted_phi  = 0.0
        attempted_psi = 0.0
        accepted_psi  = 0.0
        for k in range(0, 10):
            for residue_i in range(0,len(molecule.residues)):
                
                if molecule.fixed_residues[residue_i] != 0:
                    # se nao houver 
                    pass
                
                else:
                    j = residue_i
                    
                    if j == 1:
                        j = 2

                    if j in molecule.fixed_residues:
                        pass

                    else:
                        #--------------------------------------#
                        #                P S I                 #
                        #--------------------------------------#

                        attempted_psi += 1

                        theta = random.uniform(-1 * angle_range, angle_range)
                        rotate_backbone(molecule = molecule, 
                                            resi = j,
                                            bond = 'PSI', 
                                           theta = theta,
                                           steps = 1)
                                           

                        energy = molecule.energy()
                        if energy <= initial_energy:
                            save_XYZ_to_file(molecule, trajectory)
                            initial_energy = energy
                            initial_coordinates = molecule.get_coordinates_from_system()
                            accepted_psi += 1

                        else:
                            DG = (energy - initial_energy)
                            Px = math.exp(-1 * DG / (Kb * temp))
                            X = random.uniform(0, 1)
                            if X <= Px:
                                save_XYZ_to_file(molecule, trajectory)
                                initial_energy = energy
                                initial_coordinates = molecule.get_coordinates_from_system()
                                accepted_psi += 1

                            else:
                                molecule.import_coordinates_to_system(
                                    initial_coordinates)
                        
                        #--------------------------------------#
                        #                P H I                 #
                        #--------------------------------------#
                        theta = random.uniform(-1 * angle_range, angle_range)
                        rotate_backbone(molecule  = molecule, 
                                            resi  = j       ,
                                            bond  = 'PHI'   , 
                                            theta = theta   , 
                                            steps = 1)
                                            
                        energy = molecule.energy()

                        attempted_phi += 1

                        if energy <= initial_energy:
                            save_XYZ_to_file(molecule, trajectory)
                            initial_energy = energy
                            initial_coordinates = molecule.get_coordinates_from_system()
                            accepted_phi += 1
                        else:
                            DG = (energy - initial_energy)
                            Px = math.exp(-1 * DG / (Kb * temp))
                            X = random.uniform(0, 1)

                            if X <= Px:
                                save_XYZ_to_file(molecule, trajectory)
                                initial_energy = energy
                                initial_coordinates = molecule.get_coordinates_from_system()
                                accepted_phi += 1
                                # print 'ACCEPTED','temp:', temp,
                                # 'energy',energy,'Px',Px,'X',X
                            else:
                                # print 'Rejected'
                                # print  initial_coordinates
                                molecule.import_coordinates_to_system(
                                    initial_coordinates)
                                # print 'REJECTED','temp:', temp, 'energy',energy,'Px',Px,'X',X
                        #n = n +1

        print 'temp: = ', temp, 'energy = ', initial_energy, 'acceptance ratio (phi) =', (accepted_phi / attempted_phi), 'acceptance ratio (psi) =', (accepted_psi / attempted_psi), 'attempted_phi = ', attempted_phi
        temp = float(initial_T) * math.exp(-1 * gamma * i)
'''
