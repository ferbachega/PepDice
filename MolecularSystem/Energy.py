from ParseNaMDLog import ParseNaMDLog
from ParseGMXLog  import ParseGMXLog
from ParseAMBERLog import ParseAMBERLog
from CRDFiles      import save_CRD_to_file
from Geometry    import  distance_ab, computePhiPsi
import math 
import os
import sys


# --------- printing data --------- 
#if log:
from pprint import pprint
#    pprint(energy_list)
## ---------------------------------


class AB_ENERGY:
    """ Class doc """
    
    def __init__ (self):
        """ Class initialiser """
        pass
                        
    def compute_atomi_atomj_AB_energy (self, atom_i, atom_j, R_ab = False):
        """ Function doc """ #-23085572255.9
        #A    = self.AB_model_energy_parameters['epsilon'][0]
        #B    = 1
        
        if atom_i.AB + atom_j.AB == 'AA':   # apolares
            #C = 1.0
            #A = 10.0
            
            A = self.AB_model_energy_parameters['epsilon'][0]
            C = self.AB_model_energy_parameters['C'][0]

        
        elif atom_i.AB + atom_j.AB == 'BB': # polares
            #C = 1.0
            #A = 5
            
            A = self.AB_model_energy_parameters['epsilon'][2]
            C = self.AB_model_energy_parameters['C'][0]
        
        else:
            #C = -0.5
            #A = 5
            
            C = self.AB_model_energy_parameters['C'][1]
            A = self.AB_model_energy_parameters['epsilon'][2]

            
        #atom_i.sigma = 3.8
        #atom_j.sigma = 3.8

        sigma_ab = (atom_i.sigma_ab * atom_j.sigma_ab)**0.5
        E_ab = A*( (sigma_ab/R_ab)**12 -  C*(sigma_ab/R_ab)**6)
        return E_ab


    def compute_AB_energy (self, cutoff = 999):
        """ Function doc """ 
        '''
        ARG =  -4.5
        LYS =  -3.9  AB_ENERGY        =     -0.0000496
        ASN =  -3.5
        ASP =  -3.5
        GLU =  -3.5
        GLN =  -3.5
        HIS =  -3.2
        PRO =  -1.6
        TYR =  -1.3
        TRP =  -0.9
        SER =  -0.8
        THR =  -0.7
        GLY =  -0.4
        ALA =   1.8
        MET =   1.9
        CYS =   2.5
        PHE =   2.8
        LEU =   3.8
        VAL =   4.2
        ILE =   4.5
        
        Kyte J, Doolittle RF (May 1982). "A simple method for displaying the hydropathic character of a protein". 
        Journal of Molecular Biology.157.
        '''

        
        
        hydropathic_table = {
                            'ARG' : -4.5 / 4.5, #-4.5,
                            'LYS' : -3.9 / 4.5, #-3.9,
                            'ASN' : -3.5 / 4.5, #-3.5,
                            'ASP' : -3.5 / 4.5, #-3.5,
                            'GLU' : -3.5 / 4.5, #-3.5,
                            'GLN' : -3.5 / 4.5, #-3.5,
                            'HIS' : -3.2 / 4.5, #-3.2,
                            'HIE' : -3.2 / 4.5, #-3.2,
                            'PRO' : -1.6 / 4.5, #-1.6,
                            'TYR' : -1.3 / 4.5, #-1.3,
                            'TRP' : -0.9 / 4.5, #-0.9,
                            'SER' : -0.8 / 4.5, #-0.8,
                            'THR' : -0.7 / 4.5, #-0.7,
                            'GLY' : -0.4 / 4.5, #-0.4,
                            'ALA' :  1.8 / 4.5, # 1.8,
                            'MET' :  1.9 / 4.5, # 1.9,
                            'CYS' :  2.5 / 4.5, # 2.5,
                            'PHE' :  2.8 / 4.5, # 2.8,
                            'LEU' :  3.8 / 4.5, # 3.8,
                            'VAL' :  4.2 / 4.5, # 4.2,
                            'ILE' :  4.5 / 4.5, # 4.5,
                            }


        hydropathic_table_AB = {
                            'ARG' : 'B', #-4.5,  
                            'LYS' : 'B', #-3.9,
                            'ASN' : 'B', #-3.5,
                            'ASP' : 'B', #-3.5,
                            'GLU' : 'B', #-3.5,
                            'GLN' : 'B', #-3.5,
                            'HIS' : 'B', #-3.2,
                            'HIE' : 'B', #-3.2,

                            'PRO' : 'B', #-1.6,
                            'TYR' : 'B', #-1.3,
                            'TRP' : 'B', #-0.9,
                            'SER' : 'B', #-0.8,
                            'THR' : 'B', #-0.7,
                            'GLY' : 'B', #-0.4,
                            'ALA' : 'A', # 1.8,
                            'MET' : 'A', # 1.9,
                            'CYS' : 'A', # 2.5,
                            'PHE' : 'A', # 2.8,
                            'LEU' : 'A', # 3.8,
                            'VAL' : 'A', # 4.2,
                            'ILE' : 'A', # 4.5,
                            }

        total_E = 0

        atom_i = None
        atom_J = None
        for index_i in range(0, len(self.residues)):
            for index_j in range(index_i+2, len(self.residues)):
                
                name_i = self.residues[index_i].name
                for atom in self.residues[index_i].atoms:
                    if atom.name == 'CA':
                        atom_i    = atom  
                        atom_i.hydropathic = hydropathic_table[name_i]
                        atom_i.AB          = hydropathic_table_AB[name_i]
                        
                        #mass_i = mass_table[name_i]
                        
                        
                
                name_j = self.residues[index_j].name
                for atom in self.residues[index_j].atoms:
                    if atom.name == 'CA':
                        atom_j = atom  
                        atom_j.hydropathic = hydropathic_table[name_j]
                        atom_j.AB          = hydropathic_table_AB[name_j]
                
                        #mass_j = mass_table[name_j]



                
                R_ab = distance_ab (atom_i, atom_j)
                E    = self.compute_atomi_atomj_AB_energy (atom_i, atom_j, R_ab)
                #print index_i, name_i, hydropathic_table[name_i], atom_i.pos , index_j, name_j, hydropathic_table[name_j], atom_j.pos, 'distance_ij: ', distance_ab (atom_i, atom_j), compute_ij_CalphaModel_vdw (atom_i, atom_j, R_ab)
                
                #print '%4i %5s %10.6f %4i %5s %10.6f %10.4f %20.15f' %(index_i, name_i, hydropathic_table[name_i],  index_j, name_j, hydropathic_table[name_j], distance_ab (atom_i, atom_j), E)
                total_E += E
                
                
        #atom_i.sigma = 3.8
        #atom_j.sigma = 3.8
        #
        #
        #for distance in range(1,500):
        #    pass
       
        return total_E


class AMBER_ENERGY:
    """ Class doc """
    
    def __init__ (self):
        """ Class initialiser """
        pass

    def compute_AMBER_energy (self, pn = 1, log = None):
        """ Function doc """
        # transformar numa funcao
        write_AMBER_input_file(molecule = self, Type='energy', pn= pn, parameters = self.amber_single_point_parammeters)
        save_CRD_to_file      (molecule = self, filename='SinglePoint'+str(pn)+'.crd')
        
        os.system('sander -O -i SinglePoint'+str(pn)+'.in -c SinglePoint'+str(pn)+'.crd -o SinglePoint'+str(pn)+'.log -p ' + self.top)
        
        
        energy_list = ParseAMBERLog('SinglePoint'+str(pn)+'.log', log=log)

       
        
        if energy_list["VDWAALS"] == None:
            energy_list["VDWAALS"] = None
            return None, energy_list
        
        else:
            pass
        
        energy = 0
        
        for energy_conponent in energy_list:
            if energy_list[energy_conponent] == None:
                return None, energy_list
            else:
                energy += energy_list[energy_conponent]
        
        return energy, energy_list

    
class RG_GIRATION:
    """ Class doc """
    
    def __init__ (self):
        """ Class initialiser """
        pass

    
    def compute_center_of_mass (self):
        mass_table= {
                        'ALA'  :  71.03711 ,
                        'ARG'  :  156.10111,
                        'ASN'  :  114.04293,
                        'ASP'  :  115.02694,
                        'CYS'  :  103.00919,
                        'GLU'  :  129.04259,
                        'GLN'  :  128.05858,
                        'GLY'  :  57.02146 ,
                        
                        'HIS'  :  137.05891,
                        'HIE'  :  137.05891,
                        'HID'  :  137.05891,
                        
                        'ILE'  :  113.08406,
                        'LEU'  :  113.08406,
                        'LYS'  :  128.09496,
                        'MET'  :  131.04049,
                        'PHE'  :  147.06841,
                        'PRO'  :  97.05276 ,
                        'SER'  :  87.03203 ,
                        'THR'  :  101.04768,
                        'TRP'  :  186.07931,
                        'TYR'  :  163.06333,
                        'VAL'  :  99.06841 }
        total_E = 0
        atom_i = None
        atom_J = None
        mass_center = [0.0, 0.0, 0.0]
        mass_sum    =  0.0
        
        for index_i in range(0, len(self.residues)):
            
            name_i = self.residues[index_i].name
            for atom in self.residues[index_i].atoms:
                if atom.name == 'CA':
                    atom_i = atom  
                    mass_i = mass_table[name_i]
                    self.residues[index_i].mass = mass_i
                    
            mass_center[0] += mass_i*atom_i.pos[0]
            mass_center[1] += mass_i*atom_i.pos[1]
            mass_center[2] += mass_i*atom_i.pos[2]
            mass_sum   += mass_i
        
        mass_center[0] = mass_center[0]/mass_sum
        mass_center[1] = mass_center[1]/mass_sum
        mass_center[2] = mass_center[2]/mass_sum
        return mass_center
    

    def compute_R_gy_Calpha (self):
        """ """
        energy = 0.0
        mass_center =  self.compute_center_of_mass()
        
        Rc = ((mass_center[0]**2+mass_center[1]**2 + mass_center[2]**2)**0.5)
        
        Rg = 0.0
        mass_sum = 0
        for index_i in range(0, len(self.residues)):
            
            name_i = self.residues[index_i].name
            
            for atom in self.residues[index_i].atoms:
                
                if atom.name == 'CA':
                
                    atom_i = atom  
            
            
            Ra        = ((atom_i.pos[0]**2 + atom_i.pos[1]**2 + atom_i.pos[2]**2)**0.5)
            
            Rg       += self.residues[index_i].mass * ((Ra - Rc)**2)
            mass_sum += self.residues[index_i].mass
            
        Rg = Rg/mass_sum
        return Rg
        
        
class CONTACT_ENERGY:
    """ Class doc """
    
    def __init__ (self):
        """ Class initialiser """
        pass
    def compute_CONTACT_energy (self, log = False, cutoff = 6.0):
        """ Function doc """
        energy = 0.0
        
        for index_i in range(0, len(self.residues)):
            for index_j in range(index_i+2, len(self.residues)):
                
                name_i = self.residues[index_i].name
                
                for atom in self.residues[index_i].atoms:
                    if atom.name == 'CA':
                        atom_i    = atom  
                name_j = self.residues[index_j].name
                
                for atom in self.residues[index_j].atoms:
                    if atom.name == 'CA':
                        atom_j = atom  
                
                R_ab = distance_ab (atom_i, atom_j)
                #print index_i, name_i, index_j, name_j, R_ab
                
                
                #if self.cmap[index_i][index_j] != 0:
                #    print index_i, name_i, index_j, name_j, 'beep'
                
                if R_ab <= cutoff:
                    #print 'R_ab <= cutoff', R_ab , cutoff
                    
                    #se houver contato
                    if self.cmap[index_i][index_j] != 0:
                        #print 'beep'
                        #print energy
                        #verifica se o contato eh valido - segundo a matrix de contato
                        
                        energy += -1
                        
                        #if log:
                        #print 'beep'
                        #print index_i, name_i, index_j, name_j, R_ab, 'contact', self.cmap[index_i][index_j]
                else:
                    #print index_i, name_i, index_j, name_j, R_ab, 'NO CONTACT contact', self.cmap[index_i][index_j]
                    pass    
        
        #if log:
        #    print 'total E:', energy
        
        return energy
        

class GEOMETRY_ENERGY:
    """ Class doc """
    
    def __init__ (self):
        """ Class initialiser """
        pass
    
    def compute_SS_energy (self, log = False):
        """ Function doc """
        
        if log:
            text = '%-14s ' %('RESIDUE')
            text+= '%-7s  ' %('PHI')  
            text+= '%-14s  ' %('Restrainted-PHI')  
            text+= '%-14s  ' %('DELTA PHI')  
            text+= '%-14s  ' %('weight')  
            text+= '%-14s  ' %('phi_energy')  

            
            text+= '%-7s  ' %('PSI')  
            text+= '%-14s  ' %('Restrainted-PSI')  
            text+= '%-14s  ' %('DELTA PSI') 
            text+= '%-14s  ' %('weight')  
            text+= '%-14s  ' %('psi_energy')  
            
            text+= '%-14s  ' %('TOTAL ENERGY')  

            print text

        
        total_energy = 0
        energy_list  = []        
        
        for i in range (0, len(self.residues)):
            
            phi_final_angle = computePhiPsi (molecule=self, resi=i, bond='PHI')
            psi_final_angle = computePhiPsi (molecule=self, resi=i, bond='PSI')
            #ome_final_angle = computePhiPsi (molecule=self, resi=i, bond='OMEGA')

            
            
            phi_energy = 0
            psi_energy = 0
            Kd = 1.0
            
            if phi_final_angle != None:
                if self.residues[i].phi_restraint_angle != None:
                    delta_phi      = phi_final_angle - (self.residues[i].phi_restraint_angle)
                    K_phi =  Kd * self.residues[i].phi_restraint_weight
                    phi_energy = K_phi*(1 - math.cos(math.radians(delta_phi)))
                    #print phi_energy
                else:
                    delta_phi = False
                    self.residues[i].phi_restraint_angle = 0
                
                
            if psi_final_angle != None:
                if self.residues[i].psi_restraint_angle != None:
                    delta_psi     =  psi_final_angle - (self.residues[i].psi_restraint_angle) 
                    K_psi         =  Kd *self.residues[i].psi_restraint_weight
                    psi_energy    = K_psi*(1 - math.cos(math.radians(delta_psi)))
                else:
                    delta_psi = False
                    self.residues[i].psi_restraint_angle = 0
            


            energy     =  phi_energy + psi_energy 
            total_energy += energy

            #if i ==1:
            #    print phi_final_angle, phi_energy, self.residues[i].phi_restraint_weight, psi_final_angle, psi_energy, self.residues[i].psi_restraint_weight,  energy, total_energy


            #if i ==1:
            #    print i, phi_final_angle, energy, psi_energy, psi_energy

            #energy_list.append([self.residues[i].name, phi_energy, psi_energy, energy])
            
            '''
            if log:

                text = '%s '    %(self.residues[i].name)
                text+= '%15.5f' %(phi_final_angle)  
                text+= '%15.5f' %(self.residues[i].phi_restraint_angle)  
                text+= '%15.5f' %(delta_phi)  
                text+= '%15.5f' %(self.residues[i].phi_restraint_weight)  
                text+= '%15.5f' %(phi_energy)  

               
                text+= '%15.5f' %(psi_final_angle)  
                text+= '%15.5f' %(self.residues[i].psi_restraint_angle)  
                text+= '%15.5f' %(delta_psi)  
                text+= '%15.5f' %(self.residues[i].psi_restraint_weight)  
                text+= '%15.5f' %(psi_energy)  
                
                text+= '%15.5f' %(energy)  
               
                print text
            '''
        if log:
            print 'TOTAL ENERGY = ', total_energy
            print energy_list
        
        #print total_energy
        return total_energy
    
    
    def compute_harmonical_restraint_energies (self, log = False):
        """ Function doc """

        energy = 0
        
        if log:
            print '------------------------------- Harmonical Restraint Energies -----------------------------------'
            print 'Id(i)     name      id(j)     name         Rij        R restraint        delta R         Energy'    

        
        for restraint in self.hamonical_potential_restraint_list:
            name_i = self.residues[restraint['resi_i']].name          
            for atom in self.residues[restraint['resi_i']].atoms:
                if atom.name == restraint['atom_name_i']:
                    atom_i    = atom  

            
            name_j = self.residues[restraint['resi_j']].name
            for atom in self.residues[restraint['resi_j']].atoms:
                if atom.name == restraint['atom_name_j']:
                    atom_j = atom  
            
            R_ab    = distance_ab (atom_i, atom_j)
            delta_R = (R_ab - restraint['distance'])
            
            energy  = (delta_R**2)*restraint['Kd']
            
            if log:
                print '%-3s       %-3s        %s        %-3s %15.7f %15.7f %15.7f %15.7f ' %(restraint['resi_i'], name_i , restraint['resi_j'], name_j , R_ab, restraint['distance'], delta_R, energy)
                
                #print restraint['resi_i'], name_i , restraint['resi_j'], name_j , R_ab, restraint['distance'], delta_R, energy
        
        if log:
            print '-------------------------------------------------------------------------------------------------'

        return energy 



class Energy(AB_ENERGY, AMBER_ENERGY, RG_GIRATION, CONTACT_ENERGY, GEOMETRY_ENERGY):
    """ Class doc """

    def __init__ (self):
        """ Class initialiser """
    
        self.mass_table= {
                        'ALA'  :  71.03711 ,
                        'ARG'  :  156.10111,
                        'ASN'  :  114.04293,
                        'ASP'  :  115.02694,
                        'CYS'  :  103.00919,
                        'GLU'  :  129.04259,
                        'GLN'  :  128.05858,
                        'GLY'  :  57.02146 ,
                        'HIS'  :  137.05891,
                        'ILE'  :  113.08406,
                        'LEU'  :  113.08406,
                        'LYS'  :  128.09496,
                        'MET'  :  131.04049,
                        'PHE'  :  147.06841,
                        'PRO'  :  97.05276 ,
                        'SER'  :  87.03203 ,
                        'THR'  :  101.04768,
                        'TRP'  :  186.07931,
                        'TYR'  :  163.06333,
                        'VAL'  :  99.06841 }



    def energy(self, 
               log                       = False, 
               pn                        = 1    ,  #process number # used in multiprocess 
               external_coordinates      = False, 
               external_coordinates_type = 'pdb',
               external_coordinates_file = None , 
               
               # - - - -  novos termos - - - - -
               AMBER                     = True , 
               return_list               = False,
               ):
                    
    

        
        #defining SIZE
        self.energy_components['SIZE'][0] = len(self.residues)
        
        # AMBER energy
        #----------------------------------------------------------------------------------------
        
        sum_of_amber_coefs = 0
        
        for key in ['ANGLE'  ,'BOND'   ,'DIHED'  ,'EEL'    ,'EELEC'  ,'EGB'  ,'ESURF','NB','VDWAALS']:
            sum_of_amber_coefs += self.energy_components[key][0]
        
        if sum_of_amber_coefs != 0:       
            energy, amber_energy_list = self.compute_AMBER_energy(pn = pn, log= log)
            
            for key in  amber_energy_list:
                self.energy_components[key][0] = amber_energy_list[key]

            if energy == None:
                return None 

        else:
            energy = 0
            amber_energy_list = {}

            
        #----------------------------------------------------------------------------------------
        # AB energy 
        #----------------------------------------------------------------------------------------
        if self.energy_components['AB_ENERGY'][1] != 0.0:
            self.energy_components['AB_ENERGY'][0] = self.compute_AB_energy ()
        #----------------------------------------------------------------------------------------

        
        # CONTACT energy  
        #----------------------------------------------------------------------------------------
        if self.energy_components['CONTACT'][1] != 0.0:
            self.energy_components['CONTACT'][0] = self.compute_CONTACT_energy(log = log, cutoff = self.contact_energy_parameters['R_cutoff'])
        
        #----------------------------------------------------------------------------------------
        
        # R_GYRATION
        #----------------------------------------------------------------------------------------         
        if self.energy_components['R_GYRATION'][1] != 0.0:
            Rg = self.compute_R_gy_Calpha()
            self.energy_components['R_GYRATION'][0]= Rg
        else:
            pass
        #----------------------------------------------------------------------------------------    

        # SS_RESTRAINT
        #----------------------------------------------------------------------------------------         
        if self.energy_components['SS_RESTRAINT'][1] != 0.0:
            self.energy_components['SS_RESTRAINT'][0] = self.compute_SS_energy()
        else:
            pass
        #----------------------------------------------------------------------------------------

        # DIST_RESTRAINT
        #----------------------------------------------------------------------------------------         
        if self.energy_components['DIST_RESTRAINT'][1] != 0.0:
            self.energy_components['DIST_RESTRAINT'][0] = self.compute_harmonical_restraint_energies()
        #----------------------------------------------------------------------------------------


        operators = {}
        
        for key in self.energy_components:
            #print key
            if '/'  in key:
                
                key2 = key.split('/')
                component1  = key2[0]
                component2  = key2[1]
                self.energy_components[key][0] = self.energy_components[component1][0] / self.energy_components[component2.upper()][0]
            
            if '*'  in key:
                
                key2 = key.split('*')
                component1  = key2[0]
                component2  = key2[1]
                self.energy_components[key][0] = self.energy_components[component1][0] ** self.energy_components[component2.upper()][0]
            
            if '**'  in key:
                key2 = key.split('**')
                component1  = key2[0]
                component2  = key2[1]
                self.energy_components[key][0] = self.energy_components[component1][0] ** self.energy_components[component2.upper()][0]
        

        for key in self.energy_components:
            self.energy_components[key][0] = self.energy_components[key][0] * self.energy_components[key][1] 
            

        energy = 0      
        for component in self.energy_components:
            energy += self.energy_components[component][0] # *self.energy_components[ component  ]


        if log:
            #print '%-20s  = %15.7f   ( %15.8f )' %(key, self.energy_components[key][0], self.energy_components[key][1])
    
            text = ''
            text += '\n'
            text += '----------------------------------------------------------------------------------\n'       
            text += '                    Summary for Model "%s" ENERGY        \n' %(self.energy_model)
            text += '----------------------------------------------------------------------------------\n' 
            
            n = 1
            for key in self.energy_components:
                
                text += '%-20s = %15.7f      ' %(key, self.energy_components[key][0])#, self.energy_components[key][1])
                n +=1            
                
                if n >= 2:
                
                    text += '\n'
                
                    n = 0

            text += '----------------------------------------------------------------------------------\n'        
            text += 'Potencial Energy (total) = %18.7f \n'%(energy)
            text += '----------------------------------------------------------------------------------\n'        

            text += '\n\n'
            
            
            print text

        
        
        
        
        if return_list:
            return self.energy_components
        else:
            return energy









def write_AMBER_input_file (molecule=None, Type='energy', pn = 1, parameters = None):
    
    if parameters == None:
        parameters = {
                    'cut'     : 999.0    , 
                    'igb'     : 1        , 
                    'saltcon' : 0.2      , 
                    'gbsa'    : 1        , 
                    'rgbmax'  : 999.00000, 
                    'surften' : 0.010
                    }
        
        
    #print parameters
    #print 'cut= %4.1f , igb= %d , saltcon= %2.1f , gbsa= %d , rgbmax = %10.5f , surften = %6.3f ,' % (parameters['cut'     ],
    #                                                                                      parameters['igb'     ],
    #                                                                                      parameters['saltcon' ],
    #                                                                                      parameters['gbsa'    ],
    #                                                                                      parameters['rgbmax'  ],
    #                                                                                      parameters['surften' ])
    
    text = """  compute single-point energy 
 &cntrl
   cut=%4.1f, igb=%d, saltcon=%2.1f, gbsa=%d, rgbmax =%10.5f, surften = %6.3f,
   ntpr=1,
   nstlim = 0, dt=0.002,
   ntt=1, tempi=300.0, temp0=300.0, tautp=2.0,
   ntx=1, irest=0, ntb=0,
 &end
eof""" %(parameters['cut'     ],
         parameters['igb'     ],
         parameters['saltcon' ],
         parameters['gbsa'    ],
         parameters['rgbmax'  ],
         parameters['surften' ])
         
    output_file = open('SinglePoint'+str(pn)+'.in', "w")
    output_file.write(text)
    output_file.close()
    
def save_PDB_to_file(molecule, filename):
    with open(filename, "w") as output_file:

        text = ''
        n = 0

        for residue_i in molecule.residues:
            for atom_i in residue_i.atoms:
                #text += ("{}\t{}\n".format(atom_i.name,"\t".join([str(round(c, 2)) for c in atom_i.pos])))
                #n = n +1

                ATOM = "ATOM"
                idx = atom_i.id
                #nter              = atom_i.name

                Aname = atom_i.name
                resn = residue_i.name
                if resn == "HIE":
                    resn = "HIS"
                
                
                chainID = 'X'
                resi = str(residue_i.id)
                x = float(atom_i.pos[0])
                y = float(atom_i.pos[1])
                z = float(atom_i.pos[2])
                occ = 1.00
                tpF = 1.00
                segID = 'P2'
                element = atom_i.name[0]

                #text +=  "ATOM     1  " + atom +   " " +resn+ "  {:4d}    {:8.3f}{:8.3f}{:8.3f}  1.00  0.00          Na+\n".format(resi, float(k), float(i), float(j))

                #         ATOM, idx,  Aname," ", resn, ' ',chnID,resi,      x,     y,     z,    occ,   tpF,        segID,element," "
                text += "{:<6s}{:5d} {:<4s}{:1s}{:>3s}{:1s}{:2s}{:>3s}    {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}      {:<4s} {:2s}  {:>2s}\n".format(ATOM,
                                                                                                                                              idx,
                                                                                                                                              Aname,
                                                                                                                                              "",
                                                                                                                                              resn,
                                                                                                                                              " ",
                                                                                                                                              chainID,
                                                                                                                                              resi,
                                                                                                                                              x,
                                                                                                                                              y,
                                                                                                                                              z,
                                                                                                                                              occ,
                                                                                                                                              tpF,
                                                                                                                                              segID,
                                                                                                                                              element,
                                                                                                                                              "")


        output_file.write(text)
        output_file.close()



'''
def write_NaMD_input_file(molecule=None, Type='energy', parameters = None, pn = 1):
    """ Function doc """
    text = ''
    if molecule.ff_type == 'charmm':

        coordinates = 'SinglePoint'+str(pn)+'.pdb'
        structure = molecule.psf
        parameters = molecule.param
        paratypecharmm = 'on'

        text += '# NAMD Config file - autogenerated by NAMDgui plugin\n'
        text += '# Author: Jan Saam,  saam@charite.de                \n'
        text += '# input                                             \n'
        text += 'coordinates             ' + coordinates + '\n'
        text += 'structure               ' + structure + '\n'
        text += 'parameters              ' + parameters + '\n'
        text += 'paratypecharmm          ' + paratypecharmm + '\n\n'

        text += '# output                                                                 \n'
        text += 'set output              tmp                                              \n'
        text += 'outputname              $output                                          \n'
        text += 'dcdfile                 ${output}.dcd                                    \n'
        text += 'xstFile                 ${output}.xst                                    \n'
        text += 'dcdfreq                 50                                               \n'
        text += 'xstFreq                 50                                               \n'
        text += '                                                                         \n'
        text += 'binaryoutput            no                                               \n'
        text += 'binaryrestart           no                                               \n'
        text += 'outputEnergies          100                                              \n'
        text += 'restartfreq             1000                                             \n'
        text += '                                                                         \n'
        text += 'fixedAtoms              off                                              \n'
        text += '                                                                         \n'
        text += '# Basic dynamics                                                         \n'
        text += 'exclude                 scaled1-4                                        \n'
        text += '1-4scaling              1                                                \n'
        text += 'COMmotion               no                                               \n'
        text += 'dielectric              1.0                                              \n'
        text += '                                                                         \n'
        text += '# Simulation space partitioning                                          \n'
        text += 'switching               on                                               \n'
        text += 'switchdist              9                                                \n'
        text += 'cutoff                  10                                               \n'
        text += 'pairlistdist            12                                               \n'
        text += '                                                                         \n'
        text += '# Multiple timestepping                                                  \n'
        text += 'firsttimestep           0                                                \n'
        text += 'timestep                1                                                \n'
        text += 'stepspercycle           20                                               \n'
        text += 'nonbondedFreq           2                                                \n'
        text += 'fullElectFrequency      4                                                \n'
        text += '                                                                         \n'
        text += '# Temperature control                                                    \n'
        text += '                                                                         \n'
        text += 'set temperature         298                                              \n'
        text += 'temperature             $temperature;  # initial temperature             \n'
        text += '                                                                         \n'
        text += '                                                                         \n'
        text += 'GBIS                    on                                               \n'
        text += 'SASA                    on                                               \n'
        text += '# Scripting                                                              \n'
        
        if Type == 'energy':
            text += 'run 0                                                                    \n'
            output_file = open('SinglePoint'+str(pn)+'.namd', "w")
            output_file.write(text)
            output_file.close()
        
        if Type == 'minimize':
            text += 'minimize 100                                                             \n'
            output_file = open('Minimize'+str(pn)+'.namd', "w")
            output_file.write(text)
            output_file.close()
'''
            
            
            
#        '''
#        if self.energy_model == 'LABIO':
#            # AMBER energy modification
#            #----------------------------------------------------------------------------------------
#            if log:
#                print 'PARAMETER          RAW          Coef.'# %( energy_list['ANGLE'  ] ,  self.energy_components['ANGLE'  ])
#                print 'ANGLE     = %14.7f %14.7f' %( energy_list['ANGLE'  ] ,  self.energy_components['ANGLE'  ])
#                print 'BOND      = %14.7f %14.7f' %( energy_list['BOND'   ] ,  self.energy_components['BOND'   ])
#                print 'DIHED     = %14.7f %14.7f' %( energy_list['DIHED'  ] ,  self.energy_components['DIHED'  ])
#                print 'EEL       = %14.7f %14.7f' %( energy_list['EEL'    ] ,  self.energy_components['EEL'    ])
#                print 'EELEC     = %14.7f %14.7f' %( energy_list['EELEC'  ] ,  self.energy_components['EELEC'  ])
#                print 'EGB       = %14.7f %14.7f' %( energy_list['EGB'    ] ,  self.energy_components['EGB'    ])
#                print 'ESURF     = %14.7f %14.7f' %( energy_list['ESURF'  ] ,  self.energy_components['ESURF'  ])        
#                print 'NB        = %14.7f %14.7f' %( energy_list['NB'     ] ,  self.energy_components['NB'     ])             
#                print 'VDWAALS   = %14.7f %14.7f' %( energy_list['VDWAALS'] ,  self.energy_components['VDWAALS'])
#            
#            energy_list['ANGLE'  ] = energy_list['ANGLE'  ] * self.energy_components['ANGLE'  ] 
#            energy_list['BOND'   ] = energy_list['BOND'   ] * self.energy_components['BOND'   ] 
#            energy_list['DIHED'  ] = energy_list['DIHED'  ] * self.energy_components['DIHED'  ] 
#            energy_list['EEL'    ] = energy_list['EEL'    ] * self.energy_components['EEL'    ] 
#            energy_list['EELEC'  ] = energy_list['EELEC'  ] * self.energy_components['EELEC'  ] 
#            energy_list['EGB'    ] = energy_list['EGB'    ] * self.energy_components['EGB'    ] 
#            energy_list['ESURF'  ] = energy_list['ESURF'  ] * self.energy_components['ESURF'  ]         
#            energy_list['NB'     ] = energy_list['NB'     ] * self.energy_components['NB'     ]              
#            energy_list['VDWAALS'] = energy_list['VDWAALS'] * self.energy_components['VDWAALS']                
#            #----------------------------------------------------------------------------------------
#            
#            
#            # AB energy modification
#            #----------------------------------------------------------------------------------------
#            if log:
#                print 'AB_ENERGY = %14.7f %14.7f' %( energy_list['AB_ENERGY'] , self.energy_components['AB'])
#
#            #print 'AB_ENERGY'  , energy_list['AB_ENERGY'] , self.energy_components['AB']
#            energy_list['AB_ENERGY'] = energy_list['AB_ENERGY'] * self.energy_components['AB']
#            #----------------------------------------------------------------------------------------
#            
#            
#            # CONTACT energy modification 
#            #----------------------------------------------------------------------------------------
#            if log:
#                print 'CONTACT   = %14.7f %14.7f' %(energy_list['CONTACT'],self.energy_components['CONTACT'])
#            energy_list['CONTACT']  = energy_list['CONTACT'] *self.energy_components['CONTACT']
#            #----------------------------------------------------------------------------------------
#
#            # sum of total pseudo energy - starts with a constant:
#            if log:
#                print 'CONSTANT  = %14.7f' %(self.energy_components['CONSTANT'])
#            energy_list['CONSTANT'] = self.energy_components['CONSTANT']
#            
#            
#            
#            
#            
#            if log:
#                print ''' 
#                #Energy ~ <SIZE> + <CONTACT> + <AB_ENERGY> + <ANGLE> + <BOND> + <DIHED>
#                #        + <EEL> + <EELEC> + <EGB> + <ESURF> + <NB> + <VDWAALS> + <CONSTANT>
#                 '''
#            #Y ~ <SIZE> + <contacts0> + <AB_ENERGY> + <ANGLE> + <BOND> + <DIHED>
#            # + <EEL> + <EELEC> + <EGB> + <ESURF> + <NB> + <VDWAALS> + <CONSTANT>
#            SIZE      = len(self.residues)*self.energy_components['SIZE']
#            CONTACT = energy_list['CONTACT']
#            AB_ENERGY = energy_list['AB_ENERGY']
#            ANGLE     = energy_list['ANGLE']
#            BOND      = energy_list['BOND']
#            DIHED     = energy_list['DIHED']
#            EEL       = energy_list['EEL']
#            EELEC     = energy_list['EELEC']
#            EGB       = energy_list['EGB']
#            ESURF     = energy_list['ESURF']
#            NB        = energy_list['NB']
#            VDWAALS   = energy_list['VDWAALS']
#            intercept = energy_list['CONSTANT']
#            
#            energy = SIZE + CONTACT + AB_ENERGY + ANGLE + BOND + DIHED + EEL + EELEC + EGB + ESURF + NB + VDWAALS + intercept
#            energy = energy**(10.0/3)
#
#
#        if self.energy_model == 'iLABIO':
#            # AMBER components modification
#            #----------------------------------------------------------------------------------------
#            energy_list['ANGLE'  ] = energy_list['ANGLE'  ] * self.energy_components['ANGLE'  ] 
#            energy_list['BOND'   ] = energy_list['BOND'   ] * self.energy_components['BOND'   ] 
#            energy_list['DIHED'  ] = energy_list['DIHED'  ] * self.energy_components['DIHED'  ] 
#            energy_list['EEL'    ] = energy_list['EEL'    ] * self.energy_components['EEL'    ] 
#            energy_list['EELEC'  ] = energy_list['EELEC'  ] * self.energy_components['EELEC'  ] 
#            energy_list['EGB'    ] = energy_list['EGB'    ] * self.energy_components['EGB'    ] 
#            energy_list['ESURF'  ] = energy_list['ESURF'  ] * self.energy_components['ESURF'  ]         
#            energy_list['NB'     ] = energy_list['NB'     ] * self.energy_components['NB'     ]              
#            energy_list['VDWAALS'] = energy_list['VDWAALS'] * self.energy_components['VDWAALS']                
#            #----------------------------------------------------------------------------------------
#            
#            # AB energy modification
#            #----------------------------------------------------------------------------------------
#            energy_list['AB_ENERGY'] = energy_list['AB_ENERGY'] * self.energy_components['AB']
#            #----------------------------------------------------------------------------------------
#            
#
#            # CONTACT energy modification 
#            #----------------------------------------------------------------------------------------
#            energy_list['CONTACT'] = energy_list['CONTACT']*self.energy_components['CONTACT']
#            #----------------------------------------------------------------------------------------
#
#
#            # sum of total pseudo energy - starts with a constant:
#            energy = self.energy_components['CONSTANT']
#            
#            # sum of the components
#            for component in energy_list:
#                #print component,  energy_list[component]
#                energy += energy_list[component]
#            
#            # final component - SIZE
#            energy += len(self.residues)*self.energy_components['SIZE']
#            energy = energy**(10.0/3)
#        '''
#
