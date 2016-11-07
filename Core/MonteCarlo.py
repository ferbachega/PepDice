from MetropolisAcceptanceTest import metropolis_acceptance_test
from Geometry                 import *
from XYZFiles                 import save_XYZ_to_file

def fragment_selection (molecule           = None  ,
                        previous_fragment  = None  ,
                        fragment_rate      = None  , 
                        random             = None  ):
    
    
    fragment_acceptance = random.uniform(0, 1)
    FRAGMENTS = False
    PHIPSI    = False
    
    #----------------------------------------------#
    #                  FRAGMENTS                   #
    #----------------------------------------------#
    # (1) se o numero for menor ou igual a chance, entao um novo fragmento eh atribuido a estrutura
    if fragment_acceptance <= fragment_rate:
             
        # (2) sorteia uma posicao do alinhamento
        resi = random.randint(0, len(molecule.fragments)-1)

        # (3) enquando nao houver fragmentos para esta posicao (ou seja, lista vazia) sortear novas posicoes
        while molecule.fragments[resi] == []:
            resi = random.randint(0, len(molecule.fragments)-1) #(-1) -> nao pegar a ultima posicao
            
        # (4) sorteia um fragmento para a posicao selecionada 
        fragment_index = random.randint(0, len(molecule.fragments[resi])-1)
        fragment       = molecule.fragments[resi][fragment_index]
        
        if fragment != previous_fragment:
            return fragment, fragment_index
        
        else:
            return False
    
    else:
        return False
  

def fragment_check ():
    """ Verifies if the fragment was used in the previous step"""
    

def insert_fragment (molecule = None, fragment = None, sidechain = False):
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
            #try:
            for bond in ['CHI1','CHI2','CHI3','CHI4','CHI5']:
                #try:
                if bond in fragment[key]:
                    try:
                        set_chi_dihedral (molecule  = molecule,
                                              resi  = key,
                                              bond  = bond,
                                              angle = fragment[key][bond])
                    except:
                        pass
                #except:
                #    print 'fail', bond
            #except KeyError as error:
            #    logger.debug(error.message)
            #    logger.debug("Target: " + "".join([
            #        three_to_one(molecule.residues[i].name)
            #            if molecule.residues[i].name != 'HIE' else 'H'
            #            for i in fragment
            #    ]))
            #    logger.debug("Fragment: " + "".join([
            #        three_to_one(a['NAME']) for a in fragment.values()
            #    ]))

def rotate_backbone_attempt (molecule  = None,
                                 resi  = None,
                                 bond  = None,
                                theta  = None,
                      previous_energy  = None,
                          temperature  = None,
                                   pn  = None,
                                   Kb  = None,
                                random = None):

    rotate_backbone(molecule = molecule,
                        resi = resi    ,
                        bond = bond    ,
                       theta = theta   )

    energy = molecule.energy(pn =pn)

    if energy:
        if metropolis_acceptance_test (energy  = energy          ,
                              previous_energy  = previous_energy ,
                                  temperature  = temperature     ,
                                           Kb  = Kb              , 
                                       random  =  random):
            return energy

        else:
            return False
    else:
        return False























def monte_carlo(molecule            = None                   ,
                random              = None                   ,
                
                initial_temperature = 1000                   ,
                final_temperature   = False                  ,
                gamma               = False                  ,
                
                
                Kb                  = 1                      , # 0.0083144621               ,
                angle_range         = 60                      ,
                fragment_rate       = 1.0                    , #between 0  and 1
                fragment_sidechain  = True                   ,
                PhiPsi_rate         = 1.0                    ,
                                    
                
                simulated_annealing = None                   , # exp, linear
                cycle_size          = 10000                  ,
                number_of_cycles    = 10                     ,
                
                log_frequence       = 10                     ,
                trajectory          = 'MonteCarlo_trajectory',
                pn                  = 1                      ,
                log                 = True                   ,):
    
    
    if random == None:
        import random as random
    
    temperature          = float(initial_temperature)
    initial_temperature  = temperature
    energy = None
    
    
    
    #-------------------------------------------------#
    #                  Thermalization                 #
    #-------------------------------------------------#
    
    Temperature_Handling = 'constast'
    
    if simulated_annealing == 'exp':
        
        print '\n\n - - - simulated annealing mode - - -'
        
        if final_temperature :
        
            print 'final temperature = ', final_temperature
            
            if gamma:
                print 'ignored gamma valor   = ', gamma

            gamma = -1*(math.log(final_temperature/initial_temperature))/(number_of_cycles*cycle_size)
            
            print 'new gamma valor   = ', gamma
            Temperature_Handling = 'exponential'

        else:
            if gamma:
                print 'Using gamma (set by user) = ',  gamma
                Temperature_Handling = 'exp decay'

            else:
                print 'Error, simulated annealing mode requires a final temperature value.' 
                simulated_annealing =  False
    
    if simulated_annealing == 'linear':
        decay_factor  =  float(temperature)/float(number_of_cycles)
        Temperature_Handling = 'linear decay'

    #-------------------------------------------------#



    
    
    fragment = None 
    
    #            LOGFILE     
    logfilename = trajectory+'.log'
    infofile    = trajectory+'.info'
    
    logfile     = open(logfilename, 'a')
    infofile    = open(infofile,    'w')
    
    trajectory  = trajectory+'.xyz'
    waste       = trajectory+'.waste.xyz'  


    #            energy
    previous_energy      = molecule.energy(pn = pn)
    previous_coordinates = molecule.get_coordinates_from_system()
    previous_fragment    = None

    logfile_counter = 1
    text = []
    


    nSteps        =  0    


    #                                  printing log
    infotext = ''
    infotext += '\n'
    infotext += '\n--------------------------------- Monte Carlo Integrator Options -----------------------------'
    infotext += '\nLog Frequency               = %14d      Maximum Iterations          = %14d'    %(log_frequence        , number_of_cycles*cycle_size )
    infotext += '\nInitial temperature         = %14.5f      Final temperature           = %14.5f'%(initial_temperature  , final_temperature           ) 
    infotext += '\nTemperature Handling        = %14s      Gamma (decay factor)        = %14s'    %(Temperature_Handling , str(gamma)                  )
    infotext += '\nCycle size                  = %14d      Number of cycles            = %14d'    %(cycle_size           , number_of_cycles            )
    infotext += '\nKboltzmann                  = %14.5f      Label                       = %14s'  %(Kb                   , molecule.name               )
    infotext += '\nPhi/Psi rate                = %14.5f      Maximum angle range         = %14.3f'%(PhiPsi_rate          , angle_range                 ) 
    infotext += '\nFragment rate               = %14.5f      Fragment side chains        = %14s  '%(fragment_rate        , fragment_sidechain          ) 
    infotext += '\n----------------------------------------------------------------------------------------------'
    infotext += '\n\n'
    infotext += '\n---------------------------------------------------------------------------------'
    infotext += '\n   cycle       TEMP            E_mean          Acc_frag     Acc_phi     Acc_phi  '
    infotext += '\n---------------------------------------------------------------------------------'

    if log:
        print infotext
    infofile.write(infotext)


    for cycle in range(0, number_of_cycles):
        sum_energy = 0 
        attempted  = 0
        
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
        
        
        
        for i in range(0, cycle_size):
            nSteps +=1 

            # parametros controla a taxa de tentativas com que novos fragmentos sao testatos
            fragment_acceptance = random.uniform(0, 1)
            FRAGMENTS = False
            PHIPSI    = False
            
            
            #----------------------------------------------#
            #                 FRAGMENTS                    #
            #----------------------------------------------#
            if fragment_rate != 0.0:
                try:
                    fragment, fragment_index  = fragment_selection (molecule          = molecule         ,
                                                                    previous_fragment = previous_fragment,
                                                                    fragment_rate     = fragment_rate    ,
                                                                    random            = random           )
                except:
                    pass
                    print ' - - - fragment_selection failed - - - '
                    print 'fragment:', fragment
                
                
            
            if fragment: 
                
                attempted_fragment += 1.0
                #------- associa o fragmento selecionado com a estrutura -------
                previous_fragment = fragment
                
                insert_fragment (molecule   = molecule,
                                 fragment   = fragment,
                                 sidechain  = fragment_sidechain)
                
                #---------------------------------------------------------------
                energy = molecule.energy(pn = pn)            
                
                if energy:

                    if metropolis_acceptance_test (energy          = energy         ,
                                                   previous_energy = previous_energy,
                                                   temperature     = temperature    , 
                                                   random          = random         ,
                                                   Kb              = Kb             ):

                        save_XYZ_to_file (molecule, trajectory)
                        previous_energy      = energy
                        previous_coordinates = molecule.get_coordinates_from_system()
                        
                        accepted_fragment += 1
                        
                        FRAGMENTS = True
                        #text.append("pn: {:<3d} step: {:5d} energy: {:<20.7f}fragment: {:<3d}\n".format(pn, i, energy , fragment_index))
                        #print "pn: {:<3d} step: {:5d} energy: {:<20.7f}fragment: {:<3d}".format(pn, i, energy , fragment_index)
                    else:
                        save_XYZ_to_file (molecule, waste)                            # salva as coordenadas no descarte
                        molecule.import_coordinates_to_system (previous_coordinates)  # Restaura as coordenadas entriores 
                        #print 'fragment: ',fragment_index, energy, len(fragment), 'failed'
                else:
                    save_XYZ_to_file (molecule, waste)                                # salva as coordenadas no descarte
                    molecule.import_coordinates_to_system (previous_coordinates)      # Restaura as coordenadas entriores 



            #----------------------------------------------#
            #               PHI/PSI sampling               #
            #----------------------------------------------#

            PhiPsi_acceptance = random.uniform(0, 1)
            
            if PhiPsi_acceptance <= PhiPsi_rate:
                resi = random.randint(0, len(molecule.residues)-1)
        
                if resi in molecule.fixed_residues:
                    pass
                
                else:
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
                                                              pn = pn              ,
                                                          random = random          ,
                                                              Kb = Kb              )

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

            
            #if side_chain == True:
            #    monte_carlo_side_chain (molecule  = molecule        ,
            #                           rotamers   = rotamers        ,
            #                           random     = random          ,
            #                           initial_T  = temperature     ,
            #                           final_T    = temperature     ,
            #                           angle_range= angle_range     ,
            #                           nSteps     = side_chain_steps,
            #                           trajectory = trajectory      )
            
            
            #        got changes in energy?
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
                    text.append("pn: {:<3d} cycle: {:5d} step: {:5d} energy: {:<20.10f}".format(pn, cycle,  nSteps, energy) + frag + phipsi+'\n')
                    
                    sum_energy        += energy
                    attempted         += 1
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
            
            
           
        
        # ---------------------- Cycle report --------------------------
        if attempted_fragment != 0.0:
            acceptance_fragment = accepted_fragment / attempted_fragment
        else:
            acceptance_fragment = 0

        if attempted != 0:
            infotext = '\n%6d  %14.7f   %16.6f      %8.4f    %8.4f    %8.4f' %(cycle, temperature, sum_energy/ attempted, acceptance_fragment, (accepted_phi / attempted_phi), (accepted_psi / attempted_psi))
            if log:
                print '%6d  %14.7f   %16.6f      %8.4f    %8.4f    %8.4f' %(cycle, temperature, sum_energy/ attempted, acceptance_fragment, (accepted_phi / attempted_phi), (accepted_psi / attempted_psi))
            
            infofile.write(infotext)
        # --------------------------------------------------------------

        
        
        #-------------------Temperature Handling------------------------
        if simulated_annealing == 'exp':
            temperature = initial_temperature * math.exp(-1 * gamma * nSteps)
        
        if simulated_annealing == 'linear':
            temperature = temperature - decay_factor
        

'''
def monte_carlo(molecule            = None                   ,
                random              = None                   ,
                temperature         = 10                     ,
                Kb                  = 1                      , # 0.0083144621               ,
                angle_range         = 1                      ,
                
                #nSteps              = 10000                  ,
                fragment_rate       = 1.0                    , #between 0  and 1
                fragment_sidechain  = True                   ,
                PhiPsi_rate         = 1.0                    ,
                                    
                simulated_annealing = None                   , # exp, linear
                cycle_size          = 10000                  ,
                number_of_cycles    = 10                     ,
                
                #side_chain          = False                  ,
                #side_chain_steps    = 1000                   ,
                #rotamers            = None                   ,
                log_frequence       = 10                     ,
                trajectory          = 'MonteCarlo_trajectory',
                pn                  = 1                      ):



    if random == None:
        import random as random

    fragment = None 
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
    
    counter = 0
    
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

    
    
    if simulated_annealing == 'exp':
        pass
    
    if simulated_annealing == 'linear':
        decay_factor  =  float(temperature)/float(nSteps)
        temperature = temperature - decay_factor
    #print decay_factor
    
    nSteps = cycle_size * number_of_cycles

    for i in range(0, nSteps):
        # parametros controla a taxa de tentativas com que novos fragmentos sao testatos
        fragment_acceptance = random.uniform(0, 1)
        FRAGMENTS = False
        PHIPSI    = False
        
        
        #----------------------------------------------#
        #                 FRAGMENTS                    #
        #----------------------------------------------#
        if fragment_rate != 0.0:
            try:
                fragment, fragment_index  = fragment_selection (molecule          = molecule         ,
                                                                previous_fragment = previous_fragment,
                                                                fragment_rate     = fragment_rate    ,
                                                                random            = random           )
            except:
                pass
                print ' - - - fragment_selection failed - - - '
                #print molecule          
                #print 'previous_fragment:'
                #print previous_fragment
                
                print 'fragment:'
                print fragment    
                #print random           
            
            
        
        if fragment: 
            
            attempted_fragment += 1.0
            #------- associa o fragmento selecionado com a estrutura -------
            previous_fragment = fragment
            
            insert_fragment (molecule   = molecule,
                             fragment   = fragment,
                             sidechain  = fragment_sidechain)
            
            #---------------------------------------------------------------
            energy = molecule.energy(pn = pn)            
            
            if energy:

                if metropolis_acceptance_test (energy          = energy         ,
                                               previous_energy = previous_energy,
                                               temperature     = temperature    , 
                                               random          = random         ,
                                               Kb              = Kb             ):

                    save_XYZ_to_file (molecule, trajectory)
                    previous_energy      = energy
                    previous_coordinates = molecule.get_coordinates_from_system()
                    accepted_fragment += 1
                    FRAGMENTS = True
                    #text.append("pn: {:<3d} step: {:5d} energy: {:<20.7f}fragment: {:<3d}\n".format(pn, i, energy , fragment_index))
                    #print "pn: {:<3d} step: {:5d} energy: {:<20.7f}fragment: {:<3d}".format(pn, i, energy , fragment_index)
                else:
                    save_XYZ_to_file (molecule, waste)                            # salva as coordenadas no descarte
                    molecule.import_coordinates_to_system (previous_coordinates)  # Restaura as coordenadas entriores 
                    #print 'fragment: ',fragment_index, energy, len(fragment), 'failed'
            else:
                save_XYZ_to_file (molecule, waste)                                # salva as coordenadas no descarte
                molecule.import_coordinates_to_system (previous_coordinates)      # Restaura as coordenadas entriores 
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
                                                          pn = pn              ,
                                                      random = random          ,
                                                          Kb = Kb              )

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

        
        #if side_chain == True:
        #    monte_carlo_side_chain (molecule  = molecule        ,
        #                           rotamers   = rotamers        ,
        #                           random     = random          ,
        #                           initial_T  = temperature     ,
        #                           final_T    = temperature     ,
        #                           angle_range= angle_range     ,
        #                           nSteps     = side_chain_steps,
        #                           trajectory = trajectory      )
        
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
        print i, 'energy: ', previous_energy, 'temperature', temperature
        
        
        
        if simulated_annealing == 'exp':
            pass
        
        if simulated_annealing == 'linear':
            temperature = temperature - decay_factor

    return {'pn':pn, 'energy': previous_energy, 'coords': previous_coordinates, 'temperature': temperature }

'''







































def fragment_attempt (molecule           = None                       ,
                      temperature        = 1000                       ,
                      Kb                 = 1                          , # 0.0083144621               ,
                      fragment_rate      = 1.0                        , #between 0  and 1
                      sidechain          = False                      ,
                      
                      
                      
                      #log_frequence      = 10                         ,
                      #trajectory         = 'MonteCarlo_trajectory.xyz',
                      pn                 = 1                          ,
                      
                      random             = None                       ):
    
   
    """ Function doc """
    fragment = fragment_selection (molecule = molecule,
                                   random   = random  )
        
    if fragment != previous_fragment:
        #------- associa o fragmento selecionado com a estrutura -------
        previous_fragment = fragment
        insert_fragment (molecule  = molecule,
                         fragment  = fragment,
                         sidechain = fragment_sidechain)
        #---------------------------------------------------------------
        energy = molecule.energy(pn = pn)
        
        if energy:
            if metropolis_acceptance_test (energy = energy         ,
                                  previous_energy = previous_energy,
                                      temperature = temperature    , 
                                      random      = random         ,
                                      Kb          = Kb             ):

                save_XYZ_to_file (molecule, trajectory)
                previous_energy      = energy
                previous_coordinates = molecule.get_coordinates_from_system()
                return True
                #accepted_fragment += 1
                #FRAGMENTS = True
                #text.append("pn: {:<3d} step: {:5d} energy: {:<20.7f}fragment: {:<3d}\n".format(pn, i, energy , fragment_index))
                #print "pn: {:<3d} step: {:5d} energy: {:<20.7f}fragment: {:<3d}".format(pn, i, energy , fragment_index)
            else:
                save_XYZ_to_file (molecule, waste)                            # salva as coordenadas no descarte
                molecule.import_coordinates_to_system (previous_coordinates)  # Restaura as coordenadas entriores 
                #print 'fragment: ',fragment_index, energy, len(fragment), 'failed'
                return False

        else:
            save_XYZ_to_file (molecule, waste)                                # salva as coordenadas no descarte
            molecule.import_coordinates_to_system (previous_coordinates)      # Restaura as coordenadas entriores 
            return True

    else:
        return False




'''
def testing ():
    """ Function doc """
    import random as random
    random.seed(10)
    
    E0 = -3.0
    
    for i in range(-0,1000, 10):
        print E0, i, metropolis_acceptance_test(energy =i,previous_energy = E0, random = random)
    
    
testing()
'''


