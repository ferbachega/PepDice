#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  example01.py
#  
#  Copyright 2016 Fernando Bachega <fernando@Fenrir>
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
import os                                            #
from pprint      import pprint                       #
from Molecule    import Molecule                     #
from Geometry    import *                            #
from MonteCarlo  import monte_carlo, insert_fragment #
from XYZFiles    import save_XYZ_to_file             #
from CRDFiles    import load_CRD_from_file           #
from AATorsions  import ROTAMER_LIST                 #
                                                     #
from RMSD import compute_RMSD                        #
                                                     #
#from Test        import *                           #
                                                     #
from GeometryOptimization import minimize            #
                                                     #
from random import randint                           #
                                                     #
from Energy import save_PDB_to_file                  #
                                                     #
#----------------------------------------------------#
from Fragments import import_fragments_from_pdb
import subprocess

#-------------------------------------------------------------------------------
PEPDICE = os.environ.get('PEPDICE')
PEPDICE_EXAMPLES = os.path.join(PEPDICE, 'Examples')
PEPDICE_PARAMETER= os.path.join(PEPDICE, 'Parameters')

PEPDICE_PDBS     = os.path.join(PEPDICE, 'PDBs')
PEPDICE_OUTPUTS  = os.path.join(PEPDICE, 'outputs')

os.path.join(PEPDICE_PDBS,)

#-------------------------------------------------------------------------------


def compute_torsions (system = None, log =False):
    """ Function doc """
    
    if log:
        print '''
----------------------------------------------------------
|                   TESTING TORSIONS                     |
----------------------------------------------------------
----------------------------------------------------------
| RESIDUE       PHI            PSI             OMEGA     |
----------------------------------------------------------
'''     
    torsions = []
    for i in range (0,len(system.residues)):
        phi_final_angle = computePhiPsi (molecule=system, resi=i, bond='PHI')
        psi_final_angle = computePhiPsi (molecule=system, resi=i, bond='PSI')
        ome_final_angle = computePhiPsi (molecule=system, resi=i, bond='OMEGA')

        if phi_final_angle == None:
            phi_final_angle = 0

        if psi_final_angle == None:
            psi_final_angle = 0
            
        if ome_final_angle == None:
            ome_final_angle = 0
        
        if log:
            print "%s  %15.5f  %15.5f  %15.5f" %(system.residues[i].name , phi_final_angle, psi_final_angle, ome_final_angle)
        
        
        torsions.append([phi_final_angle, 
                         psi_final_angle, 
                         ome_final_angle])
        

    if log:
        print '''
----------------------------------------------------------
\n\n''' 
    return torsions



def refold (system = None, phi_and_psi_table = None, fileout = None):
    """ Function doc """
    
    #------------------------------------------- refolding ------------------------------------------------
    for i in range (0,len(system.residues)):
        phi_final_angle = set_phi_psi_dihedral( molecule=system, resi=i, bond='PHI'  ,angle = phi_and_psi_table[i][0]) #angle = -57 )
        psi_final_angle = set_phi_psi_dihedral( molecule=system, resi=i, bond='PSI'  ,angle = phi_and_psi_table[i][1]) #angle = -47 )
        ome_final_angle = set_phi_psi_dihedral( molecule=system, resi=i, bond='OMEGA',angle = phi_and_psi_table[i][2]) #angle = -47 )
        print system.residues[i].name ,  phi_final_angle, psi_final_angle, ome_final_angle
     #   save_XYZ_to_file(system, 'refolding.xyz')

    save_PDB_to_file(system,  fileout)
    #------------------------------------------------------------------------------------------------------


def get_sequence_from_pdb (pdbcode = None):
    """ Function doc """
    pdb     = pdbcode+'.pdb'
    pdbtext =  open(pdb, 'r')
    
    sequence      = []
    sequence_code = ''
    previous_index = None
    
    
    aa_dic = { 
         'ALA':'A',
         'ARG':'R',
         'ASN':'N',
         'ASP':'D',
         'CYS':'C',
         'GLU':'E',
         'GLN':'Q',
         'GLY':'G',
         'HIS':'H',
         'ILE':'I',
         'LEU':'L',
         'LYS':'K',
         'MET':'M',
         'PHE':'F',
         'PRO':'P',
         'SER':'S',
         'THR':'T',
         'TRP':'W',
         'TYR':'Y',
         'VAL':'V'
         }
    
    
    for line in pdbtext:
        line2 = line.split()
        #print line
        if line2[0] == 'ATOM':
            #print line
            
            try:
                index = int(line2[5])
                
                if index == previous_index:
                    pass
                else:
                    sequence.append(line2[3])
                    previous_index = index
                    sequence_code += aa_dic[line2[3]]
            except:
                print line
            #print index
            #sequence[index] = line2[3]
            
        if line2[0] == 'TER' or line2[0] == 'END' :
            break
    #print sequence
    print pdbcode, sequence_code
    return sequence_code



def build_PDBcodes_dic (dic = None):
    """ Function doc """
    PDBcodes = dic
    #wget_pdbs (PDBcodes = PDBcodes)
    for pdbcode in PDBcodes:
        sequence_code  = get_sequence_from_pdb(pdbcode = pdbcode)
        PDBcodes[pdbcode] = sequence_code
    
        #print 'failed:', pdbcode

    pprint(PDBcodes)
    return PDBcodes



def read_sequence_from_seqFile(seqfile):
    sequence = open(seqfile, 'r')
    sequence = sequence.readline()
    sequence = sequence.split()
    sequence = sequence[1]
    return sequence
    
    
  

def amber_topology_angle_force_change (filein =  None, fileout = None, force = 'E+04'):
    """ Function doc """
    filein =  open(filein, 'r')
    filein2 = filein.readlines()

    text = []
    
    for line in filein2:
        if '%FLAG ANGLE_FORCE_CONSTANT                                                      ' in line:
            inicio = filein2.index(line)

        if '%FLAG ANGLE_EQUIL_VALUE                                                         ' in line:
            final = filein2.index(line)

    #print inicio, final
    for line in range(inicio, final):
        print filein2[line]
        
        if 'E+01' in filein2[line]:
            filein2[line] = filein2[line].replace('E+01', force)
            print filein2[line]


    fileout =  open(fileout, 'w')
    fileout.writelines(filein2)
    fileout.close()
    
  


def generate_extended_structures( pdbcode = None, 
                             sequence = None, 
                           phenix_opt = True, 
                            amber_opt = True,
                          pdynamo_opt = True):
    """ Function doc """
    #----------------------------------------------------------------------
    # Criando o a estrutura de interesse na forma estendida
    #----------------------------------------------------------------------
    
    system      =  Molecule() 
    system.name =  pdbcode+'_estendida'
    system.build_peptide_from_sequence (sequence    = sequence          ,
                                        _type       = 'amber'           ,
                                        force_field = 'ff03ua.labio'    ,
                                        overwrite   = True              ,
                                        )
    
    
    #----------------------------------------------------------------------
    '''Carregando o pdb com os nomes dos atomos corrigidos'''
    #----------------------------------------------------------------------
    filename = pdbcode+'_estendida'
    save_PDB_to_file          (system,    filename+'.pdb')
    system.load_PDB_to_system (filename = filename+'.pdb')
    energy_estendida  = system.energy( log = True)
    #----------------------------------------------------------------------
            
    
    #amber_topology_angle_force_change (filein =  filename+'.top', fileout = filename+'_angles.top', force = 'E+04')
    system.import_AMBER_parameters    (top       =  filename+'.top',   
                                       torsions  = os.path.join(PEPDICE_PARAMETER, 'amber/AMBER_rotamers.dat') )
    
    
    if amber_opt:
        try:
            minimize(molecule = system,
                     imin  = 1        ,
                     maxcyc= 1000     ,
                     ncyc  = 100      ,
                     cut   = 10       ,
                     rgbmax= 999      ,
                     igb   = 1        ,
                     ntb   = 0        ,
                     ntpr  = 100      ,
                     ntr   = 0        )
                     #restraintmask = ':1-50 & @CA,N,C,O=', 
                     #restraint_wt  =  50.0 
            save_PDB_to_file(system, pdbcode+'_estendida_amber_opt.pdb')
        except:
            print 'amber opt failed'
    return system


PDBcodes =  {'1GAB': ['alpha'     ],
             #'1BX4': ['alpha_beta'],
             '1UAO': ['beta'      ],
             '1LE1': ['beta'      ],
             '1CSK': ['beta'      ],
             '2GB1': ['beta'      ],
             '1E0Q': ['beta'      ],
             '1L2Y': ['alpha'     ],
             '2JOF': ['alpha'     ],
             '1RIJ': ['alpha'     ],
             '1E0N': ['beta'      ],
             '1E0L': ['beta'      ],
             '1I6C': ['beta'      ],
             '1FME': ['alpha_beta'],
             '1PSV': ['alpha_beta'],
             '1UBQ': ['alpha_beta'],
             '1WY3': ['alpha'     ],
             '1YRF': ['alpha'     ],
             '2F4K': ['alpha'     ],
             '1VII': ['alpha'     ],
             '1EI0': ['alpha'     ],
             '1ERY': ['alpha'     ],
             '2HBA': ['alpha_beta'],
             '2HEP': ['alpha'     ],
             '1RES': ['alpha'     ],
             '1BDD': ['alpha'     ],
             '1E0G': ['alpha_beta'],
             '1BDD': ['alpha'     ],
             '1DV0': ['alpha'     ],
             '1PRB': ['alpha'     ],
             '2WXC': ['alpha'     ],
            }



for pdbcode in PDBcodes:
    
    path     = os.getcwd()
    path     = os.path.join(path,'PDBs',PDBcodes[pdbcode][0], pdbcode)
    sequence = read_sequence_from_seqFile(os.path.join(path,pdbcode+'_A.seq'))
    print sequence
#try:
    print sequence
    print pdbcode
    system = generate_extended_structures(pdbcode = pdbcode, 
                                         sequence = sequence, 
                                       phenix_opt = True, 
                                        amber_opt = True,
                                      pdynamo_opt = True)

    
    reference = Molecule()
    pdbin  = os.path.join(PEPDICE_EXAMPLES,'PDBs',PDBcodes[pdbcode][0], pdbcode, pdbcode +'_A_AMBER_minimized.pdb')
    print pdbin
    reference.load_PDB_to_system (filename = pdbin)
    
    #'/home/farminf/Programas/pepdice/Examples/PDBs/alpha/1GAB/1GAB_A_AMBER.pdb')      

    torsions = compute_torsions (system = reference, log =False)
    #print torsions
    
    for angle in torsions:
        print angle[0], angle[1], angle[2]
        #angle[2]  = 180
    refold (system = system, phi_and_psi_table = torsions, fileout = pdbcode+'_A_AMBER_refold.pdb')
    
    
    minimize(molecule = system,
                         imin  = 1        ,
                         maxcyc= 1000     ,
                         ncyc  = 1000     ,
                         cut   = 10       ,
                         rgbmax= 999      ,
                         igb   = 1        ,
                         ntb   = 0        ,
                         ntpr  = 100      ,
                         ntr   = 0        )
                         #restraintmask = ':1-50 & @CA,N,C,O=', 
                         #restraint_wt  =  50.0 
    save_PDB_to_file(system, pdbcode+'_A_AMBER_refold_minimized.pdb')
#except:
#    print 'failed geo_opt_refnament ', pdbcode 
    
    






def geo_opt_refnament (PDBcodes, phenix_opt = True, amber_opt = True):
    """ Function doc """
    for pdbcode in PDBcodes:
        
        path     = os.getcwd()
        path     = os.path.join(path,'PDBs',PDBcodes[pdbcode][0], pdbcode)
        sequence = read_sequence_from_seqFile(os.path.join(path,pdbcode+'_A.seq'))
        
        try:
            print sequence
            print pdbcode
            system = generate_extended_chains(pdbcode = pdbcode, 
                                             sequence = sequence, 
                                           phenix_opt = True, 
                                            amber_opt = True,
                                          pdynamo_opt = True)
        except:
            print 'failed geo_opt_refnament ', pdbcode
            
            
        
        lista = {}
        
        for tag in ['_A.pdb','_A_AMBER_opt.pdb','_A_AMBER_pDynamoMinimization.pdb']:#,'_A_minimized.pdb']:
            
            # Importando as informacoes torcionais da estrutura nativa do pdb
            #-------------------------------------------------------------------------------
            reference = Molecule() 
            reference.load_PDB_to_system (os.path.join(path,pdbcode+tag))  
            phi_and_psi_table  = compute_torsions(system = reference, log = False)
            #-------------------------------------------------------------------------------
            lista[pdbcode+tag] = phi_and_psi_table
            print pdbcode+tag
        print lista 
        
        #refold (system            = system, 
        #        phi_and_psi_table = raw_phi_and_psi_table, 
        #        fileout           = pdbcode+'_01_refolded_raw.pdb')
        

        
        '''
        # Criando o a estrutura de interesse na forma estendida
        #----------------------------------------------------------------------
        system = Molecule() 
        system.name =  pdbcode+'_estendida'
        system.build_peptide_from_sequence (sequence    = sequence          ,
                                            _type       = 'amber'           ,
                                            force_field = 'ff03ua'          ,
                                            overwrite   = True              ,
                                            )
        #----------------------------------------------------------------------
        #Carregando o pdb com os nomes dos atomos corrigidos
        #----------------------------------------------------------------------
        filename = pdbcode+'_estendida'
        save_PDB_to_file          (system,    filename+'.pdb')
        system.load_PDB_to_system (filename = filename+'.pdb')
        energy_estendida  = system.energy( log = True)
        #----------------------------------------------------------------------
                
        
        
        if amber_opt:
            minimize(molecule = system,
                     imin  = 1        ,
                     maxcyc= 1000     ,
                     ncyc  = 100      ,
                     cut   = 10       ,
                     rgbmax= 999      ,
                     igb   = 1        ,
                     ntb   = 0        ,
                     ntpr  = 100      ,
                     ntr   = 0        )
                     #restraintmask = ':1-50 & @CA,N,C,O=', 
                     #restraint_wt  =  50.0 
            save_PDB_to_file(system, pdbcode+'_estendida_amber_opt.pdb')
        
        
        
        if phenix_opt:
            #----------------------------------------------------------------------
            # minimizacao da cadeia estendida
            pdbou = minimize_phenix(pdbin = pdbcode+'_estendida.pdb', geo = True, geofile = 'geo.in' )
            system.load_PDB_to_system (filename = pdbcode+'_estendida_minimized.pdb')
            energy_estendida_minimized  = system.energy( log = True)
            #----------------------------------------------------------------------
        '''
        

        ## Importando as informacoes torcionais da estrutura nativa do pdb
        ##-------------------------------------------------------------------------------
        #reference = Molecule() 
        #reference.load_PDB_to_system (pdbcode+'.pdb')  
        #raw_phi_and_psi_table  = compute_torsions(system = reference, log = False)
        ##-------------------------------------------------------------------------------
        #
        #
        ## Importando as informacoes torcionais da estrutura nativa minimizada pela phenix.geometry_minimization
        ##-------------------------------------------------------------------------------
        #pdbou     = minimize_phenix(pdbin = pdbcode+'.pdb', geo = True, geofile = 'geo.in' )
        #reference = Molecule() 
        #reference.load_PDB_to_system (pdbcode+'_minimized.pdb')  
        #min_phi_and_psi_table  = compute_torsions(system = reference, log = False)
        ##-------------------------------------------------------------------------------    
        #
        #
        #
        #        
        ## Atribuicao dos angulos torcionais para a cadeira estendida nao otimizada
        #refold (system = system, phi_and_psi_table = raw_phi_and_psi_table, fileout = pdbcode+'_01_refolded_raw.pdb')
        #refold (system = system, phi_and_psi_table = min_phi_and_psi_table, fileout = pdbcode+'_01_refolded_min.pdb')
        #
        #
        #
        #
        #print energy_estendida, energy_estendida_minimized
        #
        ## Atribuicao dos angulos torcionais para a cadeira estendida otimizada
        #refold (system = system, phi_and_psi_table = raw_phi_and_psi_table, fileout = pdbcode+'_02_refolded_raw.pdb')
        #refold (system = system, phi_and_psi_table = min_phi_and_psi_table, fileout = pdbcode+'_02_refolded_min.pdb')

        



















#from pBabel           import AmberCrdFile_ToCoordinates3, AmberTopologyFile_ToSystem
#from pCore            import logFile, LogFileActive, Pickle, TestCase, TestDataSet, TestReal, TextLogFileWriter, Unpickle
#from pMolecule        import NBModelFull, SystemGeometryObjectiveFunction
#from pMoleculeScripts import HardSphereIonMobilities
#molecule              = AmberTopologyFile_ToSystem  ( "1GAB_extendida.amber11.top" )
#molecule.coordinates3 = AmberCrdFile_ToCoordinates3 ( "1GAB_extendida.crd"         )

#energy1 =  system.energy( log = True)
#
#os.system('phenix.geometry_minimization teste_novos_decoys.pdb geo.in' )
#system.load_PDB_to_system ('teste_novos_decoys_minimized.pdb') 
#energy2 =  system.energy( log = True)




#---------------------------------------------------------------------------------
#import_SS_from_string (molecule = system, ss ='CCHHHHHHHHHHHHHHHHHHHHHCCCCCHHHHHHHHHHHHHHHHHHHHHHHCHHHHHHHHHHHCC')
#save_PDB_to_file(system, 'teste_novos_decoys_SS.pdb')
#---------------------------------------------------------------------------------


'''
#---------------------------------------------------------------------------------
pdb     = 'native.pdb'      
system2 = Molecule()                                                                                                                 #
system2.load_PDB_to_system(filename = pdb)

fragment = import_fragments_from_pdb (molecule  = system2, 
                                      residues  = range(0,64), 
                                      mainchain = True, 
                                      sidechain = False)
#---------------------------------------------------------------------------------





#---------------------------------------------------------------------------------
#minimize(molecule = system,
#               imin  = 1          ,
#               maxcyc= 1000        ,
#               ncyc  = 100        ,
#               cut   = 10         ,
#               rgbmax= 999        ,
#               igb   = 1          ,
#               ntb   = 0          ,
#               ntpr  = 100        ,
#               ntr   = 0          )
#               #restraintmask = ':1-50 & @CA,N,C,O=', 
#               #restraint_wt  =  50.0 
#
#---------------------------------------------------------------------------------

insert_fragment (molecule = system,
                            fragment   = fragment,
                            sidechain  = True)
system.energy(log = True)

save_PDB_to_file(system,os.path.join(PEPDICE_EXAMPLES, 'teste_novos_decoys_fragment.pdb'))

print system.energy( log = True)
'''


'''
minimize(molecule = system,
         imin  = 1          ,
         maxcyc= 1000        ,
         ncyc  = 100        ,
         cut   = 10         ,
         rgbmax= 999        ,
         igb   = 1          ,
         ntb   = 0          ,
         ntpr  = 100        ,
         ntr   = 0          )
         #restraintmask = ':1-50 & @CA,N,C,O=', 
         #restraint_wt  =  50.0 
save_PDB_to_file(system, '1GAB_extendida_MINIMIZADA.pdb')
'''



















