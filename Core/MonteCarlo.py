from MetropolisAcceptanceTest import metropolis_acceptance_test



def fragment_attempt (molecule           = None                       ,
                      temperature        = 1000                       ,
                      Kb                 = 1                          , # 0.0083144621               ,
                      fragment_rate      = 1.0                        , #between 0  and 1
                      fragment_sidechain = False                      ,
                      log_frequence      = 10                         ,
                      trajectory         = 'MonteCarlo_trajectory.xyz',
                      pn                 = 1                          ,
                      random             = None                       ):
    """ Function doc """
    
    fragment_acceptance = random.uniform(0, 1)
    FRAGMENTS = False
    PHIPSI    = False
    
    #----------------------------------------------#
    #                  FRAGMENTS                   #
    #----------------------------------------------#
    # se o numero for menor ou igual a chance, entao um novo fragmento eh atribuido a estrutura
    if fragment_acceptance <= fragment_rate:
        
        #attempted_fragment += 1
        #sorteia uma posicao do alinhamento
        resi = random.randint(0, len(molecule.fragments)-1)

        # enquando nao houver fragmentos para esta posicao (ou seja, lista vazia) sortear novas posicoes
        while molecule.fragments[resi] == []:
            resi = random.randint(0, len(molecule.fragments)-1) #(-1) -> nao pegar aultima posicao
            
        
        # sorteia um fragmento para a posicao selecionada anteriormente 
        fragment_index = random.randint(0, len(molecule.fragments[resi])-1)
        fragment       = molecule.fragments[resi][fragment_index]
        
        #print resi, len(molecule.fragments[resi]), fragment_index, molecule.fragments[resi][fragment_index].keys()


        if fragment != previous_fragment:

            #------- associa o fragmento selecionado com a estrutura -------
            previous_fragment = fragment
            insert_fragment (molecule   = molecule,
                             fragment   = fragment,
                             sidechain  = fragment_sidechain)
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







def testing ():
    """ Function doc """
    import random as random
    random.seed(10)
    
    E0 = -3.0
    
    for i in range(-0,1000, 10):
        print E0, i, metropolis_acceptance_test(energy =i,previous_energy = E0, random = random)
    
    
testing()



