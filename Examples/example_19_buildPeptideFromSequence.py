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
    #os.system('')

def minimize_phenix (pdbin = None, geo = False, geofile = None ):
    """ Function doc """
    
    if geo:
        subprocess.call(['phenix.geometry_minimization', pdbin, geofile])
    else:
        subprocess.call(['phenix.geometry_minimization', pdbin])
    
    return pdbin+'_minimized.pdb'






PDBcodes =  {#'1GAB': []    ,
             '1UAO': []    ,
             #'1LE1':[]     ,
             #'1CSK':[]     ,
             #'2GB1':[]     ,
             #'1E0Q':[]     ,
             #'1L2Y':[]     ,
             #'2JOF':[]     ,
             #'1RIJ':[]     ,
             #'1E0N':[]     ,
             #'1E0L':[]     ,
             #'1I6C':[]     ,
             #'1FME':[]     ,
             #'1PSV':[]     ,
             #'1UBQ':[]     ,
             #'1WY3':[]     ,
             #'1YRF':[]     ,
             #'2F21':[]     ,
             #'2F4K':[]     ,
             #'1VII':[]     ,
             #'1EI0':[]     ,
             #'2F21':[]     ,
             #'1ERY':[]     ,
             #'2HBA':[]     ,
             #'2HEP':[]     ,
             #'1RES':[]     ,
             #'1BDD':[]     ,
             #'1E0G':[]     ,
             #'1BDD':[]     ,
             #'1DV0':[]     ,
             #'1PRB':[]     ,
             #'2WXC':[]     ,
            }


def wget_pdb (pdbcode = None, filesInFolder = None):
    """ Function doc """
    
    if filesInFolder == None:
        filesInFolder = os.listdir('.')
    
    if pdbcode+'.pdb' in filesInFolder:
        pass
    else:
        os.system('wget https://files.rcsb.org/download/'+pdbcode+'.pdb')    
    
    
def wget_pdbs (PDBcodes = None):
    """ Function doc """
    filesInFolder = os.listdir('.')
    
    for pdbcode in PDBcodes:
        wget_pdb (pdbcode       = pdbcode      , 
                  filesInFolder = filesInFolder)



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
        try:
        
            sequence_code  = get_sequence_from_pdb(pdbcode = pdbcode)
            PDBcodes[pdbcode] = sequence_code
        
        except:
            print 'failed:', pdbcode

    pprint(PDBcodes)
    return PDBcodes


def geo_opt_refnament (PDBcodes):
    """ Function doc """
    for pdbcode in PDBcodes:
        
        
        try:
            
            # Importando as informacoes torcionais da estrutura nativa do pdb
            #-------------------------------------------------------------------------------
            reference = Molecule() 
            reference.load_PDB_to_system (pdbcode+'.pdb')  
            raw_phi_and_psi_table  = compute_torsions(system = reference, log = False)
            #-------------------------------------------------------------------------------
            
            
            # Importando as informacoes torcionais da estrutura nativa minimizada pela phenix.geometry_minimization
            #-------------------------------------------------------------------------------
            pdbou     = minimize_phenix(pdbin = pdbcode+'.pdb', geo = True, geofile = 'geo.in' )
            reference = Molecule() 
            reference.load_PDB_to_system (pdbcode+'_minimized.pdb')  
            min_phi_and_psi_table  = compute_torsions(system = reference, log = False)
            #-------------------------------------------------------------------------------    
            
            
            # Criando o a estrutura de interesse na forma estendida
            #----------------------------------------------------------------------
            system = Molecule() 
            system.name =  pdbcode+'_estendida'
            system.build_peptide_from_sequence (sequence    = PDBcodes[pdbcode] ,
                                                _type       = 'amber'           ,
                                                force_field = 'ff03ua'          ,
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
            
                    
            # Atribuicao dos angulos torcionais para a cadeira estendida nao otimizada
            refold (system = system, phi_and_psi_table = raw_phi_and_psi_table, fileout = pdbcode+'_01_refolded_raw.pdb')
            refold (system = system, phi_and_psi_table = min_phi_and_psi_table, fileout = pdbcode+'_01_refolded_min.pdb')
            
            
            #----------------------------------------------------------------------
            # minimizacao da cadeia estendida
            pdbou = minimize_phenix(pdbin = pdbcode+'_estendida.pdb', geo = True, geofile = 'geo.in' )
            system.load_PDB_to_system (filename = pdbcode+'_estendida_minimized.pdb')
            energy_estendida_minimized  = system.energy( log = True)
            #----------------------------------------------------------------------
            
            print energy_estendida, energy_estendida_minimized
            
            # Atribuicao dos angulos torcionais para a cadeira estendida otimizada
            refold (system = system, phi_and_psi_table = raw_phi_and_psi_table, fileout = pdbcode+'_02_refolded_raw.pdb')
            refold (system = system, phi_and_psi_table = min_phi_and_psi_table, fileout = pdbcode+'_02_refolded_min.pdb')
        except:
            print 'failed: ', pdbcode
        



PDBcodes = build_PDBcodes_dic (dic = PDBcodes)
geo_opt_refnament (PDBcodes)



















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



















