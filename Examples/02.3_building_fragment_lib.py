#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  SideChainRefiner_Multiprocessing.py
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

#---------------------------------------------------------------------------------------------------#
import os                                                                                           #
from pprint      import pprint                                                                      #
from Molecule    import Molecule                                                                    #
from Geometry    import *                                                                           #
from MonteCarlo  import monte_carlo, insert_fragment, insert_fragment_from_dic#,  monte_carlo_dic, MC_replica_exchange, run_MC_replica_exchange #
from XYZFiles    import save_XYZ_to_file                                                            #
from CRDFiles    import load_CRD_from_file                                                          #
from AATorsions  import ROTAMER_LIST                                                                #
                                                                                                    #
from RMSD import compute_RMSD                                                                       #
                                                                                                    #
from Test        import *                                                                           #
                                                                                                    #
from GeometryOptimization import minimize                                                           #
                                                                                                    #
from random import randint                                                                          #
                                                                                                    #
from Energy import save_PDB_to_file                                                                 #
                                                                                                    #
#---------------------------------------------------------------------------------------------------#
from SideChainRefine import optimize_side_chain
#----------------------------------------------------#
from Fragments import build_fragment_library_from_pdbs



#-------------------------------------------------------------------------------
PEPDICE = os.environ.get('PEPDICE')
PEPDICE_EXAMPLES = os.path.join(PEPDICE, 'Examples')
PEPDICE_PARAMETER= os.path.join(PEPDICE, 'Parameters')
#-------------------------------------------------------------------------------




# building a contact map - required for Contact model calculations
#-------------------------------------------------------------------------------
#from CMAP import CMAP
#cmap = CMAP(pdb = os.path.join(PEPDICE_EXAMPLES , 'LABIO_set/1I6C/1I6C_A_AMBER_minimized.pdb'), cutoff = 6.5, log = True)
#-------------------------------------------------------------------------------


PDB      = '1GAB'
sequence = 'TIDQWLLKNAKEDAIAELKKAGITSDFYFNAINKAKTVEEVNALKNEILKAHA'

# creating a new system 
system = Molecule()
system.name = PDB+'-from_sequence' 
system.build_peptide_from_sequence (
                                     sequence    = sequence,
                                     _type       = 'amber'       ,
                                     force_field = 'ff03ua.labio',
                                     overwrite   = True          ,
                                     )
system.set_energy_model('RAW')

'''
for i in range (0,len(system.residues)):
    phi_final_angle = set_phi_psi_dihedral( molecule=system, resi=i, bond='PHI'  ,angle   = -57.00) #angle = -57 )
    psi_final_angle = set_phi_psi_dihedral( molecule=system, resi=i, bond='PSI'  ,angle   = -47.00) #angle = -47 )
    #ome_final_angle = set_phi_psi_dihedral( molecule=system, resi=i, bond='OMEGA',angle   =  179.99) #angle = -47 )
     
    
    #fragment[i] = {'PHI'  : phi_and_psi_table[i][0],#phi_final_angle, 
    #               'PSI'  : phi_and_psi_table[i][1],#psi_final_angle,
    #               'OMEGA': phi_and_psi_table[i][2],#ome_final_angle}
    #               }
                   
                   
    #fragment[i]['PSI']   =psi_final_angle
    #fragment[i]['OMEGA'] =ome_final_angle
    print system.residues[i].name ,  phi_final_angle, psi_final_angle#, ome_final_angle
'''




minimize(molecule = system,
               imin  = 1          ,
               maxcyc= 1000       ,
               ncyc  = 100        ,
               cut   = 10         ,
               rgbmax= 999        ,
               igb   = 1          ,
               ntb   = 0          ,
               ntpr  = 100        ,
               ntr   = 0          )
               #restraintmask = ':1-50 & @CA,N,C,O=', 
               #restraint_wt  =  50.0               )
save_PDB_to_file(system,  os.path.join(PEPDICE_EXAMPLES , 'new_method_example02.1_fragments_extended_min.pdb'))


system.load_PDB_to_system       (filename = os.path.join( PEPDICE_EXAMPLES , 'LABIO_set/'+PDB+'/'+PDB+'_A_AMBER_minimized.pdb'     ))

minimize(molecule = system,
               imin  = 1          ,
               maxcyc= 1000       ,
               ncyc  = 100        ,
               cut   = 10         ,
               rgbmax= 999        ,
               igb   = 1          ,
               ntb   = 0          ,
               ntpr  = 100        ,
               ntr   = 0          )
               #restraintmask = ':1-50 & @CA,N,C,O=', 
               #restraint_wt  =  50.0               )
save_PDB_to_file(system,  os.path.join(PEPDICE_EXAMPLES , 'new_method_example02.1_fragments_template_min.pdb'))





torsions = []
log = True
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




system.load_PDB_to_system       (filename = os.path.join(PEPDICE_EXAMPLES , 'new_method_example02.1_fragments_extended_min.pdb'))
save_PDB_to_file(system,  os.path.join(PEPDICE_EXAMPLES , 'new_method_example02.1_fragments_before_refoled_min.pdb'))

phi_and_psi_table = torsions
fragment = {}

for i in range (0,len(system.residues)):
    phi_final_angle = set_phi_psi_dihedral( molecule=system, resi=i, bond='PHI',angle   = phi_and_psi_table[i][0]) #angle = -57 )
    psi_final_angle = set_phi_psi_dihedral( molecule=system, resi=i, bond='PSI',angle   = phi_and_psi_table[i][1]) #angle = -47 )
    ome_final_angle = set_phi_psi_dihedral( molecule=system, resi=i, bond='OMEGA',angle = phi_and_psi_table[i][2]) #angle = -47 )
    
    
    fragment[i] = {'PHI'  : phi_and_psi_table[i][0],#phi_final_angle, 
                   'PSI'  : phi_and_psi_table[i][1],#psi_final_angle,
                   'OMEGA': phi_and_psi_table[i][2],#ome_final_angle}
                   }
                   
                   
    #fragment[i]['PSI']   =psi_final_angle
    #fragment[i]['OMEGA'] =ome_final_angle
    print system.residues[i].name ,  phi_final_angle, psi_final_angle, ome_final_angle


save_PDB_to_file(system,  os.path.join(PEPDICE_EXAMPLES , 'new_method_example02.1_fragments_after_refoled_min.pdb'))

#------------------------------------------------------------------------------------------------------
print fragment



def insert_fragment_from_dic (molecule = None, fragment = {}, sidechain = False):
    """ Function doc """
    for i in fragment:
        phi_final_angle = set_phi_psi_dihedral( molecule=molecule, resi=i, bond='PHI',angle   = fragment[i]['PHI']  )
        psi_final_angle = set_phi_psi_dihedral( molecule=molecule, resi=i, bond='PSI',angle   = fragment[i]['PSI']  )
        ome_final_angle = set_phi_psi_dihedral( molecule=molecule, resi=i, bond='OMEGA',angle = fragment[i]['OMEGA'])
        print i, phi_final_angle, psi_final_angle, ome_final_angle

    for i in fragment:
        if sidechain:
            for chi in molecule.torsions[molecule.residues[i].name]:
                
                try:
                    angle = fragment[i][chi]
                    set_chi_dihedral (molecule=molecule, resi=i, bond=chi, angle = angle, log = False)
                    #fragment[i][chi] = computeCHI (molecule= molecule, resi=i, bond=chi)
                    #print i, residue.name, chi,  fragment[i][chi] 
                except:
                    pass
                     
                    #fragment[i][chi] = None
                    #print i, residue.name, chi,  fragment[i][chi] 
            
        

system.load_PDB_to_system       (filename = os.path.join(PEPDICE_EXAMPLES , 'new_method_example02.1_fragments_extended_min.pdb'))
insert_fragment_from_dic (molecule = system, fragment = fragment)
#insert_fragment (molecule = system, fragment = fragment, sidechain = True)
save_PDB_to_file(system,  os.path.join(PEPDICE_EXAMPLES , 'new_method_example02.1_fragments_after_insertion_EXT_list.pdb'))




system.load_PDB_to_system       (filename = os.path.join(PEPDICE_EXAMPLES , 'new_method_example02.1_fragments_extended_min.pdb'))
#insert_fragment_from_dic (molecule = system, fragment = fragment)
insert_fragment (molecule = system, fragment = fragment, sidechain = False)
save_PDB_to_file(system,  os.path.join(PEPDICE_EXAMPLES , 'new_method_example02.1_fragments_after_insertion_EXT_list_2.pdb'))










pdbs          = [
                    os.path.join( PEPDICE_EXAMPLES , 'LABIO_set/'+PDB+'/'+PDB+'_A_AMBER_minimized.pdb'     ),
                    #os.path.join( PEPDICE_EXAMPLES , 'LABIO_set/1I6C/decoy0_1_A_AMBER_minimized.pdb' ),
                    #os.path.join( PEPDICE_EXAMPLES , 'LABIO_set/1I6C/decoy1_1_A_AMBER_minimized.pdb' ),
                    #os.path.join( PEPDICE_EXAMPLES , 'LABIO_set/1I6C/decoy1_2_A_AMBER_minimized.pdb' ),
                    #os.path.join( PEPDICE_EXAMPLES , 'LABIO_set/1I6C/decoy1_3_A_AMBER_minimized.pdb' ),
                    #os.path.join( PEPDICE_EXAMPLES , 'LABIO_set/1I6C/decoy1_4_A_AMBER_minimized.pdb' ),
                    #os.path.join( PEPDICE_EXAMPLES , 'LABIO_set/1I6C/decoy1_5_A_AMBER_minimized.pdb' ),
                    #os.path.join( PEPDICE_EXAMPLES , 'LABIO_set/1I6C/decoy1_6_A_AMBER_minimized.pdb' ),
                    #os.path.join( PEPDICE_EXAMPLES , 'LABIO_set/1I6C/decoy1_7_A_AMBER_minimized.pdb' ),
                    #os.path.join( PEPDICE_EXAMPLES , 'LABIO_set/1I6C/decoy1_8_A_AMBER_minimized.pdb' ),
                    #os.path.join( PEPDICE_EXAMPLES , 'LABIO_set/1I6C/decoy1_9_A_AMBER_minimized.pdb' ),
                    #os.path.join( PEPDICE_EXAMPLES , 'LABIO_set/1I6C/decoy1_10_A_AMBER_minimized.pdb'),
                    #os.path.join( PEPDICE_EXAMPLES , 'LABIO_set/1I6C/decoy1_11_A_AMBER_minimized.pdb'),
                    #os.path.join( PEPDICE_EXAMPLES , 'LABIO_set/1I6C/decoy1_12_A_AMBER_minimized.pdb'),
                    #os.path.join( PEPDICE_EXAMPLES , 'LABIO_set/1I6C/decoy1_13_A_AMBER_minimized.pdb'),
                    #os.path.join( PEPDICE_EXAMPLES , 'LABIO_set/1I6C/decoy1_14_A_AMBER_minimized.pdb'),
                    #os.path.join( PEPDICE_EXAMPLES , 'LABIO_set/1I6C/decoy1_15_A_AMBER_minimized.pdb'),
                    #os.path.join( PEPDICE_EXAMPLES , 'LABIO_set/1I6C/decoy1_16_A_AMBER_minimized.pdb'),
                    #os.path.join( PEPDICE_EXAMPLES , 'LABIO_set/1I6C/decoy1_17_A_AMBER_minimized.pdb'),
                    #os.path.join( PEPDICE_EXAMPLES , 'LABIO_set/1I6C/decoy1_23_A_AMBER_minimized.pdb'),
                    #os.path.join( PEPDICE_EXAMPLES , 'LABIO_set/1I6C/decoy1_24_A_AMBER_minimized.pdb'),
                    #os.path.join( PEPDICE_EXAMPLES , 'LABIO_set/1I6C/decoy1_30_A_AMBER_minimized.pdb'),
                    #os.path.join( PEPDICE_EXAMPLES , 'LABIO_set/1I6C/decoy1_31_A_AMBER_minimized.pdb'),
                    #os.path.join( PEPDICE_EXAMPLES , 'LABIO_set/1I6C/decoy2_14_A_AMBER_minimized.pdb'),
                    #os.path.join( PEPDICE_EXAMPLES , 'LABIO_set/1I6C/decoy2_35_A_AMBER_minimized.pdb')
                    ]





template = Molecule()
template.load_PDB_to_system       (filename = os.path.join(PEPDICE_EXAMPLES , 'new_method_example02.1_fragments_template_min.pdb'))
template.import_AMBER_parameters  (top      = None, #os.path.join(PEPDICE_EXAMPLES , 'LABIO_set/1I6C/1I6C_A_AMBER.top')   ,   
                                   torsions = os.path.join(PEPDICE_PARAMETER, 'amber/AMBER_rotamers.dat') )
print template.torsions



from Fragments import *
fragment = import_fragments_from_pdb (molecule = template, residues = range(0,len(system.residues)), mainchain =True, sidechain = True)


#extended
system.load_PDB_to_system       (filename = os.path.join(PEPDICE_EXAMPLES , 'new_method_example02.1_fragments_extended_min.pdb'))
insert_fragment_from_dic (molecule = system, fragment = fragment, sidechain = True)
save_PDB_to_file(system,  os.path.join(PEPDICE_EXAMPLES , 'new_method_example02.1_fragments_after_insertion_EXT.pdb'))

##template
#system.load_PDB_to_system       (filename = os.path.join(PEPDICE_EXAMPLES , 'new_method_example02.1_fragments_after_insertion_EXT.pdb'))
#
##system.load_PDB_to_system       (filename = os.path.join(PEPDICE_EXAMPLES , 'example02.1_fragments_template_min.pdb'))
#insert_fragment (molecule = system, fragment = fragment, sidechain = True)
#save_PDB_to_file(system,  os.path.join(PEPDICE_EXAMPLES , 'new_method_example02.1_fragments_after_insertion_TEMPLATE_2.pdb'))










#print fragment

'''
fragments  = build_fragment_library_from_pdbs (
                                                molecule             = system ,
                                                frag_size            = 5      ,
                                                number_of_fragments  = 5      ,
                                                pdblist              = pdbs   ,
                                                )
'''
#pprint (fragments)







'''
n = 3
for i in range(0,3):
    
    system  = build_fragment_library_from_pdbs (
                                                molecule             = system ,
                                                frag_size            = n      ,
                                                number_of_fragments  = 30     ,
                                                pdblist              = pdbs   ,
                                                )


    import pickle

    pickle.dump(system.fragments, open( "1gab_fragments"+str(n)+".p", "wb" ) )
    n += 2
'''
