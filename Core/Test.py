import os
from pprint      import pprint
from Molecule    import Molecule
from Geometry    import *
from MonteCarlo  import monte_carlo
from XYZFiles    import save_XYZ_to_file
from CRDFiles    import load_CRD_from_file
from AATorsions  import ROTAMER_LIST


from GeometryOptimization import minimize

from random import randint

from Energy import save_PDB_to_file




#---------------BACKBONE-------------------#
def test_computeEnergy (system = None):
    print '''
----------------------------------------------------------
|                   TESTING ENERGY                       |
----------------------------------------------------------
'''                                 
    system.energy(log = True)
    print '''
----------------------------------------------------------
\n\n''' 

def test_computeTorsions (system = None, log = False):
    """ Function doc """
    lista = []
    
    if log:
        print '''
----------------------------------------------------------
|                   TESTING TORSIONS                     |
----------------------------------------------------------
----------------------------------------------------------
| RESIDUE       PHI            PSI             OMEGA     |
----------------------------------------------------------
'''     
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
            
        lista.append([phi_final_angle, psi_final_angle, ome_final_angle])
    
    if log:
        print '''
----------------------------------------------------------
\n\n
    '''
    return lista
    
    
    
def test_rotetaPhiPsiOMEGA (system = None, theta =  0.017453292519943295,  computeTorsion = False):
    """ Function doc """
    backup  = system 
    for j in range(0,10):
        system = backup
        for i in range (0,len(system.residues)):
            
            
            try:
                rotate_backbone(molecule=system, resi=i, bond='PHI'  , theta= theta*(randint(j*-1,j)))
            except:
                print 'impossible to rotate PHI'
            
            try:
                rotate_backbone(molecule=system, resi=i, bond='PSI'  , theta=theta*(randint(j*-1,j)))
            except:
                print 'impossible to rotate PSI'

            try:
                rotate_backbone(molecule=system, resi=i, bond='OMEGA', theta=theta*(randint(j*-1,j)))
            except:
                print 'impossible to rotate OMEGA'
                
                
        save_PDB_to_file(system,  '1GAB_test_rotetaPhiPsiOMEGA_'+str(j)+'_rand.pdb')
        test_computeTorsions (system = system)

def test_refolding (system = None, angle_table = None ):
    """ Function doc """
    for i in range (0,len(system.residues)):
        phi_final_angle = set_phi_psi_dihedral( molecule=system, resi=i, bond='PHI'  ,angle = angle_table[i][0]) #angle = -57 )
        psi_final_angle = set_phi_psi_dihedral( molecule=system, resi=i, bond='PSI'  ,angle = angle_table[i][1]) #angle = -47 )
        ome_final_angle = set_phi_psi_dihedral( molecule=system, resi=i, bond='OMEGA',angle = angle_table[i][2]) #angle = -47 )
        print system.residues[i].name ,  phi_final_angle, psi_final_angle, ome_final_angle
        save_XYZ_to_file(system,TRAJECTORY)
        



#---------------SIDECHAIN-------------------#
def test_computeCHI (system = None):
    """ Function doc """
    for i in range(0, len(system.residues)):
        name = system.residues[i].name
        res  = system.residues[i]
        
        for key in system.torsions[name]:
            if key in ['PHI','PSI','OMEGA']:
                print res.id, name, key, computePhiPsi(molecule=system, resi=i, bond=key)
            
            else:
                print res.id, name, key, computeCHI(molecule=system, resi=i, bond=key)
        
def test_rotate_sidechain (system =None):
    """ Function doc """
    theta =  0.17453292519943295
    for i in range(0,len(system.residues)):
        for chi in ["CHI1","CHI2","CHI3","CHI4","CHI5"]:
            for j in range (0,50,5):
                rotate_side_chain(molecule=system, resi=i, bond=chi, theta = theta)
                
                save_XYZ_to_file(system, TRAJECTORY)
            
            print i, system.residues[i].name, chi, computeCHI(molecule=system, resi=i, bond=chi)
            
def test_set_rotamers (system = None, energy = True, trajectory = False):
    """ Function doc """
    for i in range(0, len(system.residues)):
        
        name         = system.residues[i].name
        
        res          = system.residues[i]
        
        if name in ['HSD', 'HSE', 'HDP', 'HIE', 'HID']:
            name = 'HIS'

        res_rotamers =  ROTAMER_LIST[name]

        for key in res_rotamers:
            set_side_chain_rotamer(molecule=system, resi=i, rotamer=res_rotamers[key])
            
            if energy:
                print res.id, name, key#, system.energy()
            
            if trajectory:
                save_XYZ_to_file(system, trajectory)
            


