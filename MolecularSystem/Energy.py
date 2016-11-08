from ParseNaMDLog import ParseNaMDLog
from ParseGMXLog  import ParseGMXLog
from ParseAMBERLog import ParseAMBERLog
from CRDFiles      import save_CRD_to_file
from Geometry    import  distance_ab
import math 
import os
import sys


# --------- printing data --------- 
#if log:
from pprint import pprint
#    pprint(energy_list)
## ---------------------------------



class Energy:
    """ Class doc """
    
    def __init__ (self):
        """ Class initialiser """
        pass

    def compute_atomi_atomj_AB_energy (self, atom_i, atom_j, R_ab = False):
        """ Function doc """ #-23085572255.9
        A    = 100000
        B    = 1
        
        if atom_i.AB + atom_j.AB == 'AA':   # apolares
            C = 1
        elif atom_i.AB + atom_j.AB == 'BB': # polares
            C = 0.5
        else:
            C = -0.5
    
        

        sigma_ab = (atom_i.sigma * atom_j.sigma)**0.5
        E_ab = A*( (sigma_ab/R_ab)**12 -  C*(sigma_ab/R_ab)**6)
        return E_ab


    def compute_AB_energy (self, cutoff = 999):
        """ Function doc """
        '''
        ARG =  -4.5
        LYS =  -3.9
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
                name_j = self.residues[index_j].name
                for atom in self.residues[index_j].atoms:
                    if atom.name == 'CA':
                        atom_j = atom  
                        atom_j.hydropathic = hydropathic_table[name_j]
                        atom_j.AB          = hydropathic_table_AB[name_j]
                
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

    
    def compute_AMBER_energy (self, pn = 1, log = None):
        """ Function doc """
        # transformar numa funcao
        write_AMBER_input_file(molecule = self, Type='energy', pn= pn, parameters = self.energy_model_parameters)
        save_CRD_to_file      (molecule = self, filename='SinglePoint'+str(pn)+'.crd')
        
        os.system('sander -O -i SinglePoint'+str(pn)+'.in -c SinglePoint'+str(pn)+'.crd -o SinglePoint'+str(pn)+'.log -p ' + self.top)
        
        
        energy_list = ParseAMBERLog('SinglePoint'+str(pn)+'.log', log=log)

       
        
        if energy_list["VDWAALS"] == None:
            return None
        else:
            pass
        
        energy = 0
        for energy_conponent in energy_list:
            energy += energy_list[energy_conponent]
        
        return energy, energy_list

    
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
                    
    
        energy_list = {'AB_ENERGY': 0.0,
                       'CONTACT'  : 0.0,
                       'ANGLE'    : 0.0,
                       'BOND'     : 0.0,
                       'DIHED'    : 0.0,
                       'EEL'      : 0.0,
                       'EELEC'    : 0.0,
                       'EGB'      : 0.0,
                       'ESURF'    : 0.0,
                       'NB'       : 0.0,
                       'VDWAALS'  : 0.0,}
        
        
        # AMBER energy
        #----------------------------------------------------------------------------------------
        energy, energy_list = self.compute_AMBER_energy(pn = pn, log= log)
        #----------------------------------------------------------------------------------------

        # AB energy 
        #----------------------------------------------------------------------------------------
        if self.energy_model_parameters['AB'] != 0.0:
            energy_list['AB_ENERGY'] = self.compute_AB_energy ()
            #energy_list['AB_ENERGY'] = energy_list['AB_ENERGY'] * self.energy_model_parameters['AB']
        #----------------------------------------------------------------------------------------

        # CONTACT energy  
        #----------------------------------------------------------------------------------------
        if self.energy_model_parameters['CONTACT'] != 0.0:
            energy_list['CONTACT'] = self.compute_CONTACT_energy(log = log, cutoff = self.energy_model_parameters['R_contact'])
        #energy_list['CONTACT'] = energy_list['CONTACT']*self.energy_model_parameters['CONTACT']
        #----------------------------------------------------------------------------------------
    

        if self.energy_model == 'FULL':

            energy = self.energy_model_parameters['CONSTANT']
            # sum of the components
            for component in energy_list:
                #print component,  energy_list[component]
                energy += energy_list[component]

        if self.energy_model == 'LABIO':
            # AMBER energy modification
            #----------------------------------------------------------------------------------------
            if log:
                print 'PARAMETER          RAW          Coef.'# %( energy_list['ANGLE'  ] ,  self.energy_model_parameters['ANGLE'  ])
                print 'ANGLE     = %14.7f %14.7f' %( energy_list['ANGLE'  ] ,  self.energy_model_parameters['ANGLE'  ])
                print 'BOND      = %14.7f %14.7f' %( energy_list['BOND'   ] ,  self.energy_model_parameters['BOND'   ])
                print 'DIHED     = %14.7f %14.7f' %( energy_list['DIHED'  ] ,  self.energy_model_parameters['DIHED'  ])
                print 'EEL       = %14.7f %14.7f' %( energy_list['EEL'    ] ,  self.energy_model_parameters['EEL'    ])
                print 'EELEC     = %14.7f %14.7f' %( energy_list['EELEC'  ] ,  self.energy_model_parameters['EELEC'  ])
                print 'EGB       = %14.7f %14.7f' %( energy_list['EGB'    ] ,  self.energy_model_parameters['EGB'    ])
                print 'ESURF     = %14.7f %14.7f' %( energy_list['ESURF'  ] ,  self.energy_model_parameters['ESURF'  ])        
                print 'NB        = %14.7f %14.7f' %( energy_list['NB'     ] ,  self.energy_model_parameters['NB'     ])             
                print 'VDWAALS   = %14.7f %14.7f' %( energy_list['VDWAALS'] ,  self.energy_model_parameters['VDWAALS'])
            
            energy_list['ANGLE'  ] = energy_list['ANGLE'  ] * self.energy_model_parameters['ANGLE'  ] 
            energy_list['BOND'   ] = energy_list['BOND'   ] * self.energy_model_parameters['BOND'   ] 
            energy_list['DIHED'  ] = energy_list['DIHED'  ] * self.energy_model_parameters['DIHED'  ] 
            energy_list['EEL'    ] = energy_list['EEL'    ] * self.energy_model_parameters['EEL'    ] 
            energy_list['EELEC'  ] = energy_list['EELEC'  ] * self.energy_model_parameters['EELEC'  ] 
            energy_list['EGB'    ] = energy_list['EGB'    ] * self.energy_model_parameters['EGB'    ] 
            energy_list['ESURF'  ] = energy_list['ESURF'  ] * self.energy_model_parameters['ESURF'  ]         
            energy_list['NB'     ] = energy_list['NB'     ] * self.energy_model_parameters['NB'     ]              
            energy_list['VDWAALS'] = energy_list['VDWAALS'] * self.energy_model_parameters['VDWAALS']                
            #----------------------------------------------------------------------------------------
            
            
            # AB energy modification
            #----------------------------------------------------------------------------------------
            if log:
                print 'AB_ENERGY = %14.7f %14.7f' %( energy_list['AB_ENERGY'] , self.energy_model_parameters['AB'])

            #print 'AB_ENERGY'  , energy_list['AB_ENERGY'] , self.energy_model_parameters['AB']
            energy_list['AB_ENERGY'] = energy_list['AB_ENERGY'] * self.energy_model_parameters['AB']
            #----------------------------------------------------------------------------------------
            
            
            # CONTACT energy modification 
            #----------------------------------------------------------------------------------------
            if log:
                print 'CONTACT   = %14.7f %14.7f' %(energy_list['CONTACT'],self.energy_model_parameters['CONTACT'])
            energy_list['CONTACT']  = energy_list['CONTACT'] *self.energy_model_parameters['CONTACT']
            #----------------------------------------------------------------------------------------

            # sum of total pseudo energy - starts with a constant:
            if log:
                print 'CONSTANT  = %14.7f' %(self.energy_model_parameters['CONSTANT'])
            energy_list['CONSTANT'] = self.energy_model_parameters['CONSTANT']
            
            if log:
                print ''' 
                Energy ~ <SIZE> + <CONTACT> + <AB_ENERGY> + <ANGLE> + <BOND> + <DIHED>
                        + <EEL> + <EELEC> + <EGB> + <ESURF> + <NB> + <VDWAALS> + <CONSTANT>
                 '''
            #Y ~ <SIZE> + <contacts0> + <AB_ENERGY> + <ANGLE> + <BOND> + <DIHED>
            # + <EEL> + <EELEC> + <EGB> + <ESURF> + <NB> + <VDWAALS> + <CONSTANT>
            SIZE      = len(self.residues)*self.energy_model_parameters['SIZE']
            CONTACT = energy_list['CONTACT']
            AB_ENERGY = energy_list['AB_ENERGY']
            ANGLE     = energy_list['ANGLE']
            BOND      = energy_list['BOND']
            DIHED     = energy_list['DIHED']
            EEL       = energy_list['EEL']
            EELEC     = energy_list['EELEC']
            EGB       = energy_list['EGB']
            ESURF     = energy_list['ESURF']
            NB        = energy_list['NB']
            VDWAALS   = energy_list['VDWAALS']
            intercept = energy_list['CONSTANT']
            
            energy = SIZE + CONTACT + AB_ENERGY + ANGLE + BOND + DIHED + EEL + EELEC + EGB + ESURF + NB + VDWAALS + intercept
            energy = energy**(10.0/3)


        if self.energy_model == 'iLABIO':
            # AMBER components modification
            #----------------------------------------------------------------------------------------
            energy_list['ANGLE'  ] = energy_list['ANGLE'  ] * self.energy_model_parameters['ANGLE'  ] 
            energy_list['BOND'   ] = energy_list['BOND'   ] * self.energy_model_parameters['BOND'   ] 
            energy_list['DIHED'  ] = energy_list['DIHED'  ] * self.energy_model_parameters['DIHED'  ] 
            energy_list['EEL'    ] = energy_list['EEL'    ] * self.energy_model_parameters['EEL'    ] 
            energy_list['EELEC'  ] = energy_list['EELEC'  ] * self.energy_model_parameters['EELEC'  ] 
            energy_list['EGB'    ] = energy_list['EGB'    ] * self.energy_model_parameters['EGB'    ] 
            energy_list['ESURF'  ] = energy_list['ESURF'  ] * self.energy_model_parameters['ESURF'  ]         
            energy_list['NB'     ] = energy_list['NB'     ] * self.energy_model_parameters['NB'     ]              
            energy_list['VDWAALS'] = energy_list['VDWAALS'] * self.energy_model_parameters['VDWAALS']                
            #----------------------------------------------------------------------------------------
            
            # AB energy modification
            #----------------------------------------------------------------------------------------
            energy_list['AB_ENERGY'] = energy_list['AB_ENERGY'] * self.energy_model_parameters['AB']
            #----------------------------------------------------------------------------------------
            

            # CONTACT energy modification 
            #----------------------------------------------------------------------------------------
            energy_list['CONTACT'] = energy_list['CONTACT']*self.energy_model_parameters['CONTACT']
            #----------------------------------------------------------------------------------------


            # sum of total pseudo energy - starts with a constant:
            energy = self.energy_model_parameters['CONSTANT']
            
            # sum of the components
            for component in energy_list:
                #print component,  energy_list[component]
                energy += energy_list[component]
            
            # final component - SIZE
            energy += len(self.residues)*self.energy_model_parameters['SIZE']
            energy = energy**(10.0/3)



        if log:
            text = '''
--------------------------------- Summary of Energy Terms --------------------------------
Potential Energy    =   %20.10f     BOND             =   %20.10f
EEL                 =   %20.10f     ANGLE            =   %20.10f
EELEC               =   %20.10f     DIHED            =   %20.10f
EGB                 =   %20.10f     ESURF            =   %20.10f
NB                  =   %20.10f     VDWAALS          =   %20.10f
------------------------------------------------------------------------------------------
AB_ENERGY           =   %20.10f
CONTACT             =   %20.10f
------------------------------------------------------------------------------------------
     
            ''' %(energy,
                  energy_list['BOND'     ],
                  energy_list['EEL'      ],
                  energy_list['ANGLE'    ],
                  energy_list['EELEC'    ],
                  energy_list['DIHED'    ],
                  energy_list['EGB'      ],
                  energy_list['ESURF'    ],
                  energy_list['NB'       ],
                  energy_list['VDWAALS'  ],
                  energy_list['AB_ENERGY'],
                  energy_list['CONTACT'  ],
                  #energy_list['SIZE'     ],
                  )
            
            print text
      
        
        if return_list:
            return energy_list
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

        #string = "%s%s%s%s%s%s%s%s%s%s%s%s%s%s"% (line1, index, A_name, resn, chain, resi, gap, x, y, z, b, oc, gap2, atom)
        # text.append(string+'\n')
        # print text
        #output_file.write(str(n)+ "\n\n")
        output_file.write(text)
        output_file.close()

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
            

