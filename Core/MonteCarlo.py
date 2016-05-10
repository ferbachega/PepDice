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
import logging

import random
random.seed(10)
import numpy as np
import math
import os
from Geometry import *
from XYZFiles import save_XYZ_to_file
from multiprocessing import Pool
from Bio.PDB.Polypeptide import three_to_one

logger = logging.getLogger(__name__)


class TextLofFileWriter():
    """ Class doc """

    def __init__ (self):
        """ Class initialiser """
        pass



'''
--------------------------------------------------------------
               METROPOLIS / EXCHANGE ACCEPTANCE
--------------------------------------------------------------
'''

def exchange_acceptance_test (Ei = 0, Ej = 0, Kb = 1, Ti = 1, Tj = 1, log = False):
    """ Function doc
    Y. Sugita, Y. Okamotor Chemical Physics Letters 314 ( 1999 ) 141–151

    B = 1/(Kb*T)                                                     (1)

    Delta = [Bn -Bm]*(Ei - Ej) = w(x -> X') / w(x' -> X)             (2)


    ----------------------------------------------------
                  Metropolis Criterion
    ----------------------------------------------------

    w(x -> X') =  ( 1 ,            for Delta <= 0,                   (3)
                  ( exp(-Delta),   for Delta >  0,                   (4)

    """


    Bi = 1/(Kb*Ti)                                                   # (1)
    Bj = 1/(Kb*Tj)

    Delta = (Bj - Bi)*(Ei - Ej)                                      # (2)
    if log:
        print '\n\n'
        print 'Bi    = ', Bi
        print 'Bj    = ', Bj
        print 'dE    = ',Bj - Bi
        print 'Delta = ',(Ei - Ej)*(Bj - Bi)


    #-------------------------------------------------#
    #             Metropolis Criterion                #
    #-------------------------------------------------#
    #-------------------------------------------------#
    if Delta <= 0:                                                   # (3)
        if log:
            print 'Delta <= 0 -->', Delta
        return True


    else:                                                            # (4)
        Px = math.exp(-1*Delta)                                      # Px = exp(-Delta)
        X  = random.uniform(0, 1)
        if log:
            print 'Delta >  0 -->', Delta, 'exp(-Delta) =', Px, 'X = ', X
            #print 'Delta = ', Delta, 'Px =', Px, X
        return X <= Px
    #-------------------------------------------------#

def metropolis_acceptance_test (energy = None        ,
                       previous_energy = None        ,
                           temperature = 273.15      ,
                                    Kb = 0.0019872041):
    """ Function doc """
    if energy < previous_energy:
        return True
    else:
        dE = (energy - previous_energy)
        #Px = math.exp(-1 * dE / (Kb * temperature))
        Px = math.exp(-1 * dE / (temperature))
        if temperature >= 800:
            print temperature,dE, Px, Kb
        X  = random.uniform(0, 1)
        return X <= Px


'''
--------------------------------------------------------------
                       GEOMETRY ATTEMPT
--------------------------------------------------------------
'''
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
        if metropolis_acceptance_test (energy = energy          ,
                              previous_energy = previous_energy ,
                                  temperature = 273.15          ,
                                           Kb = 0.0019872041    ):
            return energy

        else:
            return False
    else:
        return False


def insert_fragment_attempt (molecule = None, fragment = None, sidechain = False):
    """ Function doc """
    if fragment is None:
        return
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
            except KeyError as error:
                logger.debug(error.message)
                logger.debug("Target: " + "".join([
                    three_to_one(molecule.residues[i].name)
                        if molecule.residues[i].name != 'HIE' else 'H'
                        for i in fragment
                ]))
                logger.debug("Fragment: " + "".join([
                    three_to_one(a['NAME']) for a in fragment.values()
                ]))



'''
--------------------------------------------------------------
                        MONTE CARLO
--------------------------------------------------------------
'''
def monte_carlo_dic (parameters):
    """
            Function doc
    Esta funcao existe para que a funcao monte_carlo, usada no metodo pool, tenha apenas
    uma variavel como entrada -> necessaria para a paralelizacao.

    """
    random.seed(parameters['pn'])
    return  monte_carlo(**parameters)


def monte_carlo(molecule           = None                       ,
                temperature        = 1000                       ,
                Kb                 = 0.0019872041               , # 0.0083144621               ,
                angle_range        = 1                          ,
                nSteps             = 10000                      ,
                fragment_rate      = 1.0                        , #between 0  and 1
                fragment_sidechain = False                      ,
                log_frequence      = 10                         ,
                PhiPsi_rate        = 1.0                        ,
                trajectory         = 'MonteCarlo_trajectory.xyz',
                pn                 = 1                          ):


    #LOGFILE    = trajectory+'.log'
    logfilename = trajectory+'.log'
    logfile     = open(logfilename, 'a')
    trajectory  = trajectory+'.xyz'
    waste       = trajectory+'.waste.xyz'  


    previous_energy      = molecule.energy(pn = pn)
    previous_coordinates = molecule.get_coordinates_from_system()
    previous_fragment    = None


    logfile_counter = 1
    text = []
    n = 0
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

        # parametros controla a taxa de tentativas com que novos fragmentos sao testatos
        fragment_acceptance = random.uniform(0, 1)



        FRAGMENTS = False
        PHIPSI    = False

        #----------------------------------------------#
        #                  FRAGMENTS                   #
        #----------------------------------------------#
        # se o numero for menor ou igual a chance, entao um novo fragmento eh atribuido a estrutura
        if fragment_acceptance <= fragment_rate:
            attempted_fragment += 1

            #sorteia uma posicao do alinhamento
            resi = random.randint(0, len(molecule.fragments)-1)

            # enquando nao haver fragmentos para esta posicao (ou seja, lista vazia) sortear novas posicoes
            while molecule.fragments[resi] == []:
                resi = random.randint(0, len(molecule.fragments)-1) #(-1) -> nao pegar aultima posicao

            fragment_index = random.randint(0, len(molecule.fragments[resi])-1)
            fragment       = molecule.fragments[resi][fragment_index]
            #print resi, len(molecule.fragments[resi]), fragment_index, molecule.fragments[resi][fragment_index].keys()


            if fragment != previous_fragment:
                previous_fragment = fragment

                insert_fragment_attempt (molecule   = molecule,
                                 fragment   = fragment,
                                 sidechain  = fragment_sidechain)

                energy = molecule.energy(pn = pn)
                if energy:
                    if metropolis_acceptance_test (energy = energy         ,
                                          previous_energy = previous_energy,
                                              temperature = temperature    ):

                        save_XYZ_to_file (molecule, trajectory)
                        previous_energy      = energy
                        previous_coordinates = molecule.get_coordinates_from_system()
                        accepted_fragment += 1
                        FRAGMENTS = True
                        #text.append("pn: {:<3d} step: {:5d} energy: {:<20.7f}fragment: {:<3d}\n".format(pn, i, energy , fragment_index))
                        #print "pn: {:<3d} step: {:5d} energy: {:<20.7f}fragment: {:<3d}".format(pn, i, energy , fragment_index)
                    else:
                        save_XYZ_to_file (molecule, waste)
                        molecule.import_coordinates_to_system (previous_coordinates)
                        #print 'fragment: ',fragment_index, energy, len(fragment), 'failed'
                else:
                    save_XYZ_to_file (molecule, waste)
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
                    theta = math.radians(theta)

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
                        PHIPSI = True
                        #text.append("pn: {:<3d} step: {:5d} energy: {:<20.7f}rotate_backbone theta: {:<6.3f}\n".format(pn, i, energy , theta*57.324))
                        #print "pn: {:<3d} step: {:5d} energy: {:<20.7f}rotate_backbone theta: {:<6.3f}".format(pn, i, energy , theta*57.324)
                    else:
                        save_XYZ_to_file (molecule, waste)
                        molecule.import_coordinates_to_system (previous_coordinates)

        if energy:
            if energy <= previous_energy:
                if FRAGMENTS:
                    frag   = "fragment: {:<4d}  ".format(fragment_index)
                else:
                    frag   = 'fragment: -     ';format('-')
                if PHIPSI:
                    phipsi = "rotate_backbone: {:<6.3f}".format(theta*57.324)
                else:
                    phipsi = 'rotate_backbone: None'
                text.append("pn: {:<3d} step: {:5d} energy: {:<20.10f}".format(pn, i, energy) + frag + phipsi+'\n')

        #---------------------------------------------#
        #                LOGFILEWRITE                 #
        #---------------------------------------------#
        logfile_counter += logfile_counter
        if logfile_counter >= log_frequence:
            logfile.writelines(text)
            logfile.close()
            text    = []
            logfile = open(logfilename, 'a')
            logfile_counter = 1
        #----------------------------------------------

    return {'pn':pn, 'energy': previous_energy, 'coords': previous_coordinates, 'temperature': temperature }


def monte_carlo_side_chain (molecule  = None,
                           initial_T  = 1000,
                           final_T    = 1   ,
                           angle_range= 20  ,
                           nSteps     = 1000,
                           trajectory='MonteCarlo_trajectory.xyz'):

    Kb              = 0.0019872041 #, #0.0083144621
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
--------------------------------------------------------------
                       REPLICA EXCHANGE
--------------------------------------------------------------
'''

def MC_replica_exchange (replicas   = [],
                         CPUs       = 1 ,
                         N_replicas = 1 ,
                         nExchanges = 1 ):


    for i in range(0, nExchanges):

        #--------------------------------------------#
        #                MONTE CARLO                 #
        #--------------------------------------------#
        p = Pool(CPUs)                               #
        results = p.map(monte_carlo_dic, replicas)   #
        #--------------------------------------------#

        #------------------------------------- Exchange -----------------------------------------#
        replica_results = {}
        for result in results:
            #print 'replica: %3i energy: %10.7f temp: %i' %(result['pn']        ,
            #                                               result['energy']    ,
            #                                               result['temperature'])#, len(result[2])
            replica_results[result['pn']] = {
                                'energy'     : result['energy'],
                                'coords'     : result['coords'],
                                'temperature': result['temperature']
                                }
        #----------------------------------------------------------------------------------------#


        print '\n\n'
        # teste feito para verificar quais sao os pares de replicas que serao trocados
        for i in range (1,len(replica_results)):
            Ei = replica_results[i]  ['energy']                                                                       #
            Ej = replica_results[i+1]['energy']                                                                       #

            Ti = replica_results[i]  ['temperature']                                                                  #
            Tj = replica_results[i+1]['temperature']                                                                  #

            COORDi = replica_results[i]  ['coords']                                                                   #
            COORDj = replica_results[i+1]['coords']                                                                   #


            Kb = 1 #0.0019872041
            Bi = 1/(Kb*Ti)                                                   # (1)
            Bj = 1/(Kb*Tj)
            Delta = (Bj - Bi)*(Ei - Ej)                                      # (2)

            # exchange acceptance test - teste de aceitacao de troca das replicas
            if exchange_acceptance_test(Ei=Ei, Ej=Ej , Kb = Kb, Ti=Ti, Tj=Tj ):
                
                # replicas comeca em 0 - pois eh uma lista
                
                replicas[i]  ['molecule'].import_coordinates_to_system (COORDi)
                replicas[i-1]['molecule'].import_coordinates_to_system (COORDj)
             
                print 'replica ',i+1, ' ----> ', i    , 'Delta = ', Delta, 'accepted'
                print 'replica ',i, ' ----> '  , i+1  , 'Delta = ', Delta, 'accepted'

            # caso contrario nao ha troca de coordenadas entre as replicar
            else:
                replicas[i-1]['molecule'].import_coordinates_to_system (COORDi)
                replicas[i]  ['molecule'].import_coordinates_to_system (COORDj)
                #print 'replica ',i,   '----> coord', i       , 'dE = ', (Ei - Ej)  '
                #print 'replica ',i+1, '----> coord', i+1     , 'dE = ', (Ei - Ej)  '
                print 'replica ',i+1, ' ----> ', i    , 'Delta = ', Delta, 'rejected'
                print 'replica ',i, ' ----> '  , i+1  , 'Delta = ', Delta, 'rejected'



        '''
        # teste feito para verificar quais sao os pares de replicas que serao trocados
        #--------------------------------------------------------------------------------------------------------------#
        test = random.randint(1,2)                                                                                     #
                                                                                                                       #
        # partindo de repĺica 1   1-2 3-4 5-6 ...                                                                      #
        if test  == 1:                                                                                                 #
            for i in range (1,len(replica_results)+1,2):                                                                      #
                if i == len(replica_results):                                                                                 #
                    pass                                                                                               #
                else:                                                                                                  #
                    Ei = replica_results[i]  ['energy']                                                                       #
                    Ej = replica_results[i+1]['energy']                                                                       #
                                                                                                                       #
                    Ti = replica_results[i]  ['temperature']                                                                  #
                    Tj = replica_results[i+1]['temperature']                                                                  #
                                                                                                                       #
                    COORDi = replica_results[i]  ['coords']                                                                   #
                    COORDj = replica_results[i+1]['coords']                                                                   #
                                                                                                                       #
                                                                                                                       #
                    # se o criteiro de troca for satisfeito:                                                           #
                    if exchange_acceptance_test(Ei=Ei, Ej=Ej , Kb = 0.0019872041, Ti=Ti, Tj=Tj ):                      #
                        replicas[i+1-1]['molecule'].import_coordinates_to_system (COORDi)                              #
                        replicas[i-1]  ['molecule'].import_coordinates_to_system (COORDj)                              #
                        print 'replica ',i+1, '----> coord', i    , 'dE = ', (Ei - Ej)                                 #
                        print 'replica ',i, '----> coord'  , i+1  , 'dE = ', (Ei - Ej)                                 #
                                                                                                                       #
                    # caso contrario nao ha troca de coordenadas entre as replicar                                     #
                    else:                                                                                              #
                        replicas[i  -1]['molecule'].import_coordinates_to_system (COORDi)                              #
                        replicas[i+1-1]['molecule'].import_coordinates_to_system (COORDj)                              #
                        print 'replica ',i,   '----> coord', i       , 'dE = ', (Ei - Ej)                              #
                        print 'replica ',i+1, '----> coord', i+1     , 'dE = ', (Ei - Ej)                              #
                                                                                                                       #
        # partindo de repĺica 2   2-3 4-5 6-7 ... 1-n                                                                  #
        else:                                                                                                          #
            for i in range (2,len(replica_results)+1,2):                                                                      #
                if i == len(replica_results):                                                                                 #
                    Ei = replica_results[1]['energy']                                                                         #
                    Ej = replica_results[i]['energy']                                                                         #
                                                                                                                       #
                    Ti = replica_results[1]['temperature']                                                                    #
                    Tj = replica_results[i]['temperature']                                                                    #
                                                                                                                       #
                    COORDi = replica_results[1]['coords']                                                                     #
                    COORDj = replica_results[i]['coords']                                                                     #
                                                                                                                       #
                    # se o criteiro de troca for satisfeito:                                                           #
                    if exchange_acceptance_test(Ei=Ei, Ej=Ej , Kb = 0.0019872041, Ti=Ti, Tj=Tj ):                      #
                        replicas[i-1]['molecule'].import_coordinates_to_system (COORDi)                                #
                        replicas[1-1]['molecule'].import_coordinates_to_system (COORDj)                                #
                        print 'replica ',i,  '----> coord',  1  , 'dE = ', (Ei - Ej)                                   #
                        print 'replica ',1, '----> coord' ,  i  , 'dE = ', (Ei - Ej)                                   #
                                                                                                                       #
                                                                                                                       #
                    # caso contrario nao ha troca de coordenadas entre as replicar                                     #
                    else:                                                                                              #
                        replicas[1-1]['molecule'].import_coordinates_to_system (COORDi)                                #
                        replicas[i-1]['molecule'].import_coordinates_to_system (COORDj)                                #
                        print 'replica ',1, '----> coord', 1   , 'dE = ', (Ei - Ej)                                    #
                        print 'replica ',i, '----> coord', i   , 'dE = ', (Ei - Ej)                                    #
                                                                                                                       #
                else:                                                                                                  #
                    Ei = replica_results[i]['energy']                                                                         #
                    Ej = replica_results[i+1]['energy']                                                                       #
                                                                                                                       #
                    Ti = replica_results[i]['temperature']                                                                    #
                    Tj = replica_results[i+1]['temperature']                                                                  #
                                                                                                                       #
                    COORDi = replica_results[i]  ['coords']                                                                   #
                    COORDj = replica_results[i+1]['coords']                                                                   #
                                                                                                                       #
                    # se o criteiro de troca for satisfeito:                                                           #
                    if exchange_acceptance_test(Ei=Ei, Ej=Ej , Kb = 0.0019872041, Ti=Ti, Tj=Tj ):                      #
                        replicas[i+1-1]['molecule'].import_coordinates_to_system (COORDi)                              #
                        replicas[i  -1]['molecule'].import_coordinates_to_system (COORDj)                              #
                        print 'replica ',i+1, '----> coord', i   , 'dE = ', (Ei - Ej)                                  #
                        print 'replica ',i  , '----> coord', i+1 , 'dE = ', (Ei - Ej)                                  #
                                                                                                                       #
                                                                                                                       #
                    # caso contrario nao ha troca de coordenadas entre as replicar                                     #
                    else:                                                                                              #
                        replicas[i  -1]['molecule'].import_coordinates_to_system (COORDi)                              #
                        replicas[i+1-1]['molecule'].import_coordinates_to_system (COORDj)                              #
                        print 'replica ',i  , '----> coord', i   , 'dE = ', (Ei - Ej)                                  #
                        print 'replica ',i+1 ,'----> coord', i+1 , 'dE = ', (Ei - Ej)                                  #
                                                                                                                       #
        #--------------------------------------------------------------------------------------------------------------#
        '''


def run_MC_replica_exchange (
                            molecule           = None         ,

                            N_replicas         = 1            , # >= number of CPUs
                            CPUs               = 1            , # Number os CPUs

                            min_temp           = 50           ,
                            max_temp           = 250          ,
                            PhiPsi_rate        = 0.1          ,
                            max_angle_range    = 5            ,
                            trajectory         = 'MC_replica_',
                            log_frequence      = 100          ,

                            Kb                 = 0.0019872041 , #0.0083144621 ,

                            nSteps             = 1000         , # (bloco de sim) numero de passos na simulacao de MC
                            nExchanges         = 5            , # numero de eventos de troca  ->  total de simulacao eh dado pelo  nExchanges x nSteps

                            fragment_rate      = 0.5          ,
                            fragment_sidechain = True         ,
                            log                = False        ,
                            overwrite          = True         ,
                            ):

    """ Function doc """



    temperature_factor = (max_temp-min_temp)/N_replicas
    #--------------------montagen das replicas-------------------------#
    replicas   = []  # lista, onde os elementos sao dicionarios        #
    for i in range(1, N_replicas + 1):                                 #
        parameters = {}                                                #
        parameters['molecule'          ] = molecule                    #
        parameters['temperature'       ] = min_temp                    #
        parameters['Kb'                ] = Kb                          #
        parameters['nSteps'            ] = nSteps                      #
        parameters['fragment_rate'     ] = fragment_rate               #
        parameters['fragment_sidechain'] = fragment_sidechain          #
        parameters['log_frequence'     ] = log_frequence               #
        parameters['PhiPsi_rate'       ] = PhiPsi_rate                 #
        parameters['angle_range'       ] = max_angle_range             #
        parameters['trajectory'        ] = trajectory +str(i)          #
        parameters['pn'                ] = i                           #
        replicas.append(parameters)                                    #
                                                                       #
        print min_temp
        min_temp += temperature_factor                                 #
        #print min_temp
                                                                       #
                                                                       #
        if overwrite:                                                  #
            try:                                                       #
                os.remove (trajectory +str(i)+'.log')                  #
            except:                                                    #
                pass                                                   #
                                                                       #
            try:                                                       #
                os.remove (trajectory +str(i)+'.xyz')                  #
            except:                                                    #
                pass                                                   #
    #------------------------------------------------------------------#

    print '----------------------------------------'
    print 'number of residues: ' , len(molecule.residues)
    print 'number of fragments:' , len(molecule.fragments)
    print 'number of replicas: ' , len(replicas)
    print 'number of CPUs:     ' , CPUs
    print 'temperature factor: ' , temperature_factor
    print '----------------------------------------'


    #monte_carlo_dic(replicas[0])

    MC_replica_exchange(replicas   = replicas  ,
                        CPUs       = CPUs      ,
                        N_replicas = N_replicas,
                        nExchanges = nExchanges)



##---------------------------------------------------------------------------------------
#for i in range(10):
#    Ei = random.uniform(-i,i)
#    Ej = random.uniform(-i,i)
#    Ti = 275
#    Tj = 300
#    print exchange_acceptance_test(Ei=Ei, Ej=Ej , Kb = 0.0019872041, Ti=Ti, Tj=Tj )
##---------------------------------------------------------------------------------------



