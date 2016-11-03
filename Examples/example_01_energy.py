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
from MonteCarlo  import monte_carlo                  #
from XYZFiles    import save_XYZ_to_file             #
from CRDFiles    import load_CRD_from_file           #
from AATorsions  import ROTAMER_LIST                 #
                                                     #
from RMSD import compute_RMSD                        #
                                                     #
from Test        import *                            #
                                                     #
from GeometryOptimization import minimize            #
                                                     #
from random import randint                           #
                                                     #
from Energy import save_PDB_to_file                  #
                                                     #
#----------------------------------------------------#
from SideChainRefine import optimize_side_chain
#----------------------------------------------------#
from Fragments import build_fragment_library_from_pdbs














#eita

#-------------------------------------------------------------------------------
PEPDICE = os.environ.get('PEPDICE')
PEPDICE_EXAMPLES = os.path.join(PEPDICE, 'Examples')
PEPDICE_PARAMETER= os.path.join(PEPDICE, 'Parameters')
#-------------------------------------------------------------------------------









#/home/fernando/programs/pepdice/Examples/PDBs/3Drobot/1I6C/1I6C_A_AMBER_minimized.pdb
#/home/fernando/programs/pepdice/Examples/PDBs/3Drobot/1I6C/1I6C_A_AMBER.top
#-------------------------------------------------------------------------------
system = Molecule() 
system.set_energy_model('amber')
system.load_PDB_to_system      (filename = os.path.join(PEPDICE_EXAMPLES , 'PDBs/3Drobot/1I6C/1I6C_A_AMBER_minimized.pdb'   )   )   
system.import_AMBER_parameters(top       = os.path.join(PEPDICE_EXAMPLES , 'PDBs/3Drobot/1I6C/1I6C_A_AMBER.top')   ,   
                                torsions = os.path.join(PEPDICE_PARAMETER, 'amber/AMBER_rotamers.dat') )   
#-------------------------------------------------------------------------------
system.AB = 1


print system.energy(log       = True,
                    AB_energy = True)

print '\n\n\n'
system.set_energy_model('Calpha')

print system.energy(log       = True,
                    AB_energy = True)










def build_contact_map_from_pdb (pdb = None, cutoff = 6.0, log= False):
    """ Function doc """
    import numpy as np
    from Geometry    import  distance_ab
    
    #pdb = os.path.join(PEPDICE_EXAMPLES , 'PDBs/alpha/1GAB/1GAB_A_AMBER_minimized.pdb'   )
    
    system = Molecule() 
    system.load_PDB_to_system      (filename = pdb   )   
    cmap = 0*np.random.rand(len(system.residues),len(system.residues))
    
    n_contacts =  0
    
    for index_i in range(0, len(system.residues)):
        for index_j in range(index_i+2, len(system.residues)):
            
            name_i = system.residues[index_i].name
            
            for atom in system.residues[index_i].atoms:
                if atom.name == 'CA':
                    atom_i    = atom  
            name_j = system.residues[index_j].name
            
            for atom in system.residues[index_j].atoms:
                if atom.name == 'CA':
                    atom_j = atom  
            R_ab = distance_ab (atom_i, atom_j)
            
            
            if R_ab <= cutoff:
                cmap[index_i][index_j] = 1
                n_contacts += 1
                if log:
                    print index_i, name_i, index_j, name_j, R_ab, 'contact'
            else:
                pass
    print '\ncmap matrix:'
    print cmap
    
    print '\n---------------------------------------------'
    print 'total number of residues = ', len(system.residues)
    print 'Cutoff size              = ', cutoff
    print 'number of contacts       = ', n_contacts
    print '---------------------------------------------\n'
    return cmap


cmap = build_contact_map_from_pdb (pdb  = os.path.join(PEPDICE_EXAMPLES , 'PDBs/3Drobot/1I6C/1I6C_A_AMBER_minimized.pdb'   ))
system.import_CMAP_matrix (cmap = cmap)



system.build_peptide_from_sequence ( sequence    = 'KLPPGWEKRMSRSSGRVYYFNHITNASQWERPSGNSSSG',
                                     _type       = 'amber'    ,
                                     force_field = 'ff03ua.labio'   ,
                                     overwrite   = True       ,
                                     )
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





system.set_energy_model('Contact')
#system.set_energy_model('LPFSF')
#system.set_energy_model('amber')
#system.set_energy_model('Calpha')

#print 'cmap energy:', 

system.energy(log =True)

folder = os.path.join(PEPDICE_EXAMPLES , 'PDBs/3Drobot/1I6C/')
decoys =  {
'1I6C_A_AMBER_minimized.pdb'    :          0.00,
'decoy0_1_A_AMBER_minimized.pdb':          0.23,
'decoy1_1_A_AMBER_minimized.pdb':          1.37,
'decoy1_2_A_AMBER_minimized.pdb':          1.41,
'decoy1_3_A_AMBER_minimized.pdb':          1.15,
'decoy1_4_A_AMBER_minimized.pdb':          1.39,
'decoy1_5_A_AMBER_minimized.pdb':          1.78,
'decoy1_6_A_AMBER_minimized.pdb':          1.60,
'decoy1_7_A_AMBER_minimized.pdb':          1.24,
'decoy1_8_A_AMBER_minimized.pdb':          1.95,
'decoy1_9_A_AMBER_minimized.pdb':          1.96,
'decoy1_10_A_AMBER_minimized.pdb':         1.40,
'decoy1_11_A_AMBER_minimized.pdb':         1.62,
'decoy1_20_A_AMBER_minimized.pdb':         2.85,
'decoy1_22_A_AMBER_minimized.pdb':         2.81,
'decoy1_19_A_AMBER_minimized.pdb':         2.40,
'decoy1_12_A_AMBER_minimized.pdb':         2.79,
'decoy1_13_A_AMBER_minimized.pdb':         1.99,
'decoy1_14_A_AMBER_minimized.pdb':         2.55,
'decoy1_15_A_AMBER_minimized.pdb':         2.93,
'decoy1_16_A_AMBER_minimized.pdb':         2.64,
'decoy1_17_A_AMBER_minimized.pdb':         2.18,
'decoy1_18_A_AMBER_minimized.pdb':         2.21,
'decoy1_21_A_AMBER_minimized.pdb':         2.54,
'decoy1_23_A_AMBER_minimized.pdb':         2.65,
'decoy1_24_A_AMBER_minimized.pdb':         2.82,
'decoy1_31_A_AMBER_minimized.pdb':         3.02,
'decoy1_30_A_AMBER_minimized.pdb':         3.39,
'decoy2_14_A_AMBER_minimized.pdb':         5.84,
'decoy2_35_A_AMBER_minimized.pdb':         6.84,
'decoy4_1_A_AMBER_minimized.pdb':          4.96,
'decoy4_8_A_AMBER_minimized.pdb':          5.90,
'decoy4_27_A_AMBER_minimized.pdb':         6.28,
'decoy4_44_A_AMBER_minimized.pdb':         6.66,
'decoy4_92_A_AMBER_minimized.pdb':         7.89,
'decoy4_156_A_AMBER_minimized.pdb':        8.11,
'decoy5_2_A_AMBER_minimized.pdb':          2.88,
'decoy5_1_A_AMBER_minimized.pdb':          2.66,
'decoy5_3_A_AMBER_minimized.pdb':          2.90,
'decoy5_7_A_AMBER_minimized.pdb':          3.48,
'decoy5_21_A_AMBER_minimized.pdb':         4.23,
'decoy5_23_A_AMBER_minimized.pdb':         5.74,
'decoy5_100_A_AMBER_minimized.pdb':        6.42,
'decoy5_120_A_AMBER_minimized.pdb':        6.98,
'decoy6_22_A_AMBER_minimized.pdb':         7.02,
'decoy6_46_A_AMBER_minimized.pdb':         7.67,
'decoy6_85_A_AMBER_minimized.pdb':         7.97,
'decoy6_144_A_AMBER_minimized.pdb':        8.98,
'decoy7_3_A_AMBER_minimized.pdb':          4.98,
'decoy7_8_A_AMBER_minimized.pdb':          5.77,
'decoy7_22_A_AMBER_minimized.pdb':         6.97,
'decoy7_21_A_AMBER_minimized.pdb':         6.26,
'decoy7_30_A_AMBER_minimized.pdb':         7.08,
'decoy7_38_A_AMBER_minimized.pdb':         8.08,
'decoy8_118_A_AMBER_minimized.pdb':       10.23,
'decoy9_1_A_AMBER_minimized.pdb':          3.68,
'decoy9_6_A_AMBER_minimized.pdb':          4.36,
'decoy9_22_A_AMBER_minimized.pdb':         5.60,
'decoy10_10_A_AMBER_minimized.pdb':        3.69,
'decoy10_13_A_AMBER_minimized.pdb':        3.96,
'decoy10_1_A_AMBER_minimized.pdb':         4.04,
'decoy10_40_A_AMBER_minimized.pdb':        3.91,
'decoy10_23_A_AMBER_minimized.pdb':        4.66,
'decoy10_24_A_AMBER_minimized.pdb':        4.20,
'decoy10_87_A_AMBER_minimized.pdb':        5.93,
'decoy10_76_A_AMBER_minimized.pdb':        5.82,
'decoy10_141_A_AMBER_minimized.pdb':       6.58,
'decoy10_262_A_AMBER_minimized.pdb':       9.14,
'decoy11_1_A_AMBER_minimized.pdb':         3.84,
'decoy11_21_A_AMBER_minimized.pdb':        4.52,
'decoy11_11_A_AMBER_minimized.pdb':        4.54,
'decoy11_50_A_AMBER_minimized.pdb':        5.31,
'decoy11_173_A_AMBER_minimized.pdb':      10.44,
'decoy12_86_A_AMBER_minimized.pdb':       10.21,
'decoy13_10_A_AMBER_minimized.pdb':        6.86,
'decoy13_41_A_AMBER_minimized.pdb':        7.96,
'decoy13_51_A_AMBER_minimized.pdb':        7.96,
'decoy13_96_A_AMBER_minimized.pdb':        8.09,
'decoy13_122_A_AMBER_minimized.pdb':       8.09,
'decoy13_143_A_AMBER_minimized.pdb':       9.04,
'decoy14_2_A_AMBER_minimized.pdb':         4.92,
'decoy14_7_A_AMBER_minimized.pdb':         5.72,
'decoy15_1_A_AMBER_minimized.pdb':         2.56,
'decoy15_2_A_AMBER_minimized.pdb':         3.90,
'decoy15_14_A_AMBER_minimized.pdb':        4.83,
'decoy15_29_A_AMBER_minimized.pdb':        5.34,
'decoy15_31_A_AMBER_minimized.pdb':        7.01,
'decoy15_75_A_AMBER_minimized.pdb':        7.91,
'decoy15_130_A_AMBER_minimized.pdb':       8.27,
'decoy16_1_A_AMBER_minimized.pdb':         2.39,
'decoy16_2_A_AMBER_minimized.pdb':         2.42,
'decoy16_3_A_AMBER_minimized.pdb':         2.43,
'decoy16_6_A_AMBER_minimized.pdb':         3.80,
'decoy16_8_A_AMBER_minimized.pdb':         3.13,
'decoy16_18_A_AMBER_minimized.pdb':        3.86,
'decoy16_4_A_AMBER_minimized.pdb':         3.76,
'decoy16_28_A_AMBER_minimized.pdb':        4.89,
'decoy17_5_A_AMBER_minimized.pdb':         5.90,
'decoy17_14_A_AMBER_minimized.pdb':        6.14,
'decoy17_25_A_AMBER_minimized.pdb':        7.76,
'decoy17_72_A_AMBER_minimized.pdb':        8.35,
'decoy17_140_A_AMBER_minimized.pdb':       9.77,
'decoy17_92_A_AMBER_minimized.pdb':        9.99,
'decoy17_160_A_AMBER_minimized.pdb':      10.24,
'decoy18_1_A_AMBER_minimized.pdb':         8.02,
'decoy18_72_A_AMBER_minimized.pdb':        8.44,
'decoy18_98_A_AMBER_minimized.pdb':        9.35,
'decoy18_160_A_AMBER_minimized.pdb':      10.43,
'decoy18_136_A_AMBER_minimized.pdb':      10.38,
'decoy18_148_A_AMBER_minimized.pdb':      10.29,
'decoy18_116_A_AMBER_minimized.pdb':      10.45,
'decoy19_6_A_AMBER_minimized.pdb':         4.84,
'decoy19_19_A_AMBER_minimized.pdb':        5.53,
'decoy21_1_A_AMBER_minimized.pdb':         3.75,
'decoy21_11_A_AMBER_minimized.pdb':        4.59,
'decoy21_24_A_AMBER_minimized.pdb':        5.89,
'decoy21_30_A_AMBER_minimized.pdb':        6.98,
'decoy21_51_A_AMBER_minimized.pdb':        7.92,
'decoy21_92_A_AMBER_minimized.pdb':        8.50,
'decoy21_141_A_AMBER_minimized.pdb':       9.16,
'decoy22_1_A_AMBER_minimized.pdb':         2.96,
'decoy22_12_A_AMBER_minimized.pdb':        3.75,
'decoy22_31_A_AMBER_minimized.pdb':        4.49,
'decoy22_34_A_AMBER_minimized.pdb':        4.76,
'decoy22_73_A_AMBER_minimized.pdb':        5.78,
'decoy22_68_A_AMBER_minimized.pdb':        5.17,
'decoy23_1_A_AMBER_minimized.pdb':         4.58,
'decoy23_7_A_AMBER_minimized.pdb':         5.36,
'decoy23_69_A_AMBER_minimized.pdb':        6.28,
'decoy23_34_A_AMBER_minimized.pdb':        6.84,
'decoy23_104_A_AMBER_minimized.pdb':       7.61,
'decoy23_163_A_AMBER_minimized.pdb':       8.33,
'decoy23_200_A_AMBER_minimized.pdb':       9.58,
'decoy24_2_A_AMBER_minimized.pdb':         2.43,
'decoy24_3_A_AMBER_minimized.pdb':         2.56,
'decoy24_5_A_AMBER_minimized.pdb':         2.15,
'decoy24_1_A_AMBER_minimized.pdb':         2.42,
'decoy24_4_A_AMBER_minimized.pdb':         2.75,
'decoy24_31_A_AMBER_minimized.pdb':        3.38,
'decoy24_11_A_AMBER_minimized.pdb':        3.63,
'decoy24_14_A_AMBER_minimized.pdb':        3.96,
'decoy24_23_A_AMBER_minimized.pdb':        3.77,
'decoy24_25_A_AMBER_minimized.pdb':        3.77,
'decoy24_63_A_AMBER_minimized.pdb':        4.36,
'decoy24_74_A_AMBER_minimized.pdb':        4.31,
'decoy24_82_A_AMBER_minimized.pdb':        4.82,
'decoy24_38_A_AMBER_minimized.pdb':        4.30,
'decoy24_78_A_AMBER_minimized.pdb':        4.74,
'decoy24_96_A_AMBER_minimized.pdb':        5.22,
'decoy24_135_A_AMBER_minimized.pdb':       5.15,
'decoy24_86_A_AMBER_minimized.pdb':        5.36,
'decoy24_138_A_AMBER_minimized.pdb':       6.11,
'decoy24_155_A_AMBER_minimized.pdb':       8.49,
'decoy25_36_A_AMBER_minimized.pdb':        3.26,
'decoy25_14_A_AMBER_minimized.pdb':        3.56,
'decoy25_17_A_AMBER_minimized.pdb':        3.49,
'decoy25_5_A_AMBER_minimized.pdb':         3.30,
'decoy25_8_A_AMBER_minimized.pdb':         3.29,
'decoy25_10_A_AMBER_minimized.pdb':        3.97,
'decoy25_11_A_AMBER_minimized.pdb':        3.81,
'decoy25_16_A_AMBER_minimized.pdb':        3.29,
'decoy25_21_A_AMBER_minimized.pdb':        3.92,
'decoy25_40_A_AMBER_minimized.pdb':        4.73,
'decoy26_5_A_AMBER_minimized.pdb':         6.81,
'decoy26_110_A_AMBER_minimized.pdb':       9.88,
'decoy26_133_A_AMBER_minimized.pdb':       9.99,
'decoy27_1_A_AMBER_minimized.pdb':         5.45,
'decoy27_137_A_AMBER_minimized.pdb':       9.99,
'decoy29_4_A_AMBER_minimized.pdb':         5.82,
'decoy29_124_A_AMBER_minimized.pdb':       9.98,
'decoy29_133_A_AMBER_minimized.pdb':      10.32,
'decoy32_3_A_AMBER_minimized.pdb':         6.87,
'decoy32_57_A_AMBER_minimized.pdb':        7.62,
'decoy32_78_A_AMBER_minimized.pdb':        8.88,
'decoy32_118_A_AMBER_minimized.pdb':       8.98,
'decoy32_162_A_AMBER_minimized.pdb':      10.19,
'decoy33_1_A_AMBER_minimized.pdb':         3.71,
'decoy33_30_A_AMBER_minimized.pdb':        4.97,
'decoy33_67_A_AMBER_minimized.pdb':        5.09,
'decoy33_78_A_AMBER_minimized.pdb':        5.10,
'decoy33_136_A_AMBER_minimized.pdb':       6.02,
'decoy33_101_A_AMBER_minimized.pdb':       7.06,
'decoy33_106_A_AMBER_minimized.pdb':       6.81,
'decoy33_147_A_AMBER_minimized.pdb':       7.18,
'decoy33_197_A_AMBER_minimized.pdb':       8.11,
'decoy33_216_A_AMBER_minimized.pdb':       9.65,
'decoy34_3_A_AMBER_minimized.pdb':         8.07,
'decoy34_14_A_AMBER_minimized.pdb':        8.89,
'decoy34_70_A_AMBER_minimized.pdb':        9.61,
'decoy34_74_A_AMBER_minimized.pdb':        9.19,
'decoy34_129_A_AMBER_minimized.pdb':      10.70,
'decoy34_128_A_AMBER_minimized.pdb':      10.25,
'decoy34_144_A_AMBER_minimized.pdb':      10.61,
'decoy36_162_A_AMBER_minimized.pdb':      11.30,
'decoy36_163_A_AMBER_minimized.pdb':      11.00,
'decoy37_5_A_AMBER_minimized.pdb':         7.82,
'decoy37_43_A_AMBER_minimized.pdb':        8.26,
'decoy37_73_A_AMBER_minimized.pdb':        9.82,
'decoy37_112_A_AMBER_minimized.pdb':      10.09,
'decoy39_37_A_AMBER_minimized.pdb':        7.80,
'decoy39_67_A_AMBER_minimized.pdb':        8.52,
'decoy39_113_A_AMBER_minimized.pdb':       9.19,
'decoy39_122_A_AMBER_minimized.pdb':      10.01,
'decoy40_1_A_AMBER_minimized.pdb':         4.64,
'decoy40_11_A_AMBER_minimized.pdb':        5.77,
'decoy40_26_A_AMBER_minimized.pdb':        6.59,
'decoy40_73_A_AMBER_minimized.pdb':        7.29,
'decoy40_137_A_AMBER_minimized.pdb':       8.40,
'decoy40_152_A_AMBER_minimized.pdb':       9.15,
'decoy41_139_A_AMBER_minimized.pdb':      10.54,
'decoy43_5_A_AMBER_minimized.pdb':         4.58,
'decoy45_64_A_AMBER_minimized.pdb':       11.52,
'decoy45_65_A_AMBER_minimized.pdb':       11.40,
'decoy46_1_A_AMBER_minimized.pdb':         4.10,
'decoy48_1_A_AMBER_minimized.pdb':         4.82,
'decoy48_5_A_AMBER_minimized.pdb':         5.17,
'decoy49_128_A_AMBER_minimized.pdb':      10.85,
'decoy49_131_A_AMBER_minimized.pdb':      11.15,
'decoy50_9_A_AMBER_minimized.pdb':         7.90,
'decoy50_55_A_AMBER_minimized.pdb':        8.94,
'decoy50_118_A_AMBER_minimized.pdb':       9.28,
'decoy50_111_A_AMBER_minimized.pdb':       9.16,
'decoy50_109_A_AMBER_minimized.pdb':       9.39,
'decoy50_151_A_AMBER_minimized.pdb':      10.03,
'decoy54_31_A_AMBER_minimized.pdb':        6.85,
'decoy54_48_A_AMBER_minimized.pdb':        7.91,
'decoy54_104_A_AMBER_minimized.pdb':       8.53,
'decoy54_173_A_AMBER_minimized.pdb':       8.95,
'decoy54_147_A_AMBER_minimized.pdb':       9.10,
'decoy54_192_A_AMBER_minimized.pdb':      10.18,
'decoy54_199_A_AMBER_minimized.pdb':      11.81,
'decoy54_202_A_AMBER_minimized.pdb':      11.90,
'decoy54_205_A_AMBER_minimized.pdb':      11.81,
'decoy54_208_A_AMBER_minimized.pdb':      11.86,
'decoy54_209_A_AMBER_minimized.pdb':      11.91,
'decoy54_211_A_AMBER_minimized.pdb':      11.94,
'decoy54_213_A_AMBER_minimized.pdb':      11.98,
'decoy54_198_A_AMBER_minimized.pdb':      11.89,
'decoy54_201_A_AMBER_minimized.pdb':      11.83,
'decoy54_203_A_AMBER_minimized.pdb':      11.98,
'decoy54_204_A_AMBER_minimized.pdb':      11.69,
'decoy54_206_A_AMBER_minimized.pdb':      11.95,
'decoy54_207_A_AMBER_minimized.pdb':      12.02,
'decoy54_200_A_AMBER_minimized.pdb':      11.20,
'decoy54_210_A_AMBER_minimized.pdb':      11.02,
'decoy55_164_A_AMBER_minimized.pdb':      10.49,
'decoy55_154_A_AMBER_minimized.pdb':      10.40,
'decoy58_3_A_AMBER_minimized.pdb':         6.69,
'decoy58_36_A_AMBER_minimized.pdb':        8.00,
'decoy58_19_A_AMBER_minimized.pdb':        7.83,
'decoy58_96_A_AMBER_minimized.pdb':        8.46,
'decoy58_58_A_AMBER_minimized.pdb':        8.28,
'decoy58_105_A_AMBER_minimized.pdb':       9.30,
'decoy59_48_A_AMBER_minimized.pdb':        6.82,
'decoy59_83_A_AMBER_minimized.pdb':        7.28,
'decoy59_91_A_AMBER_minimized.pdb':        7.87,
'decoy59_105_A_AMBER_minimized.pdb':       8.68,
'decoy59_138_A_AMBER_minimized.pdb':       8.58,
'decoy59_157_A_AMBER_minimized.pdb':       9.36,
'decoy59_166_A_AMBER_minimized.pdb':      10.03,
'decoy60_151_A_AMBER_minimized.pdb':      10.38,
'decoy61_6_A_AMBER_minimized.pdb':         7.07,
'decoy61_8_A_AMBER_minimized.pdb':         6.94,
'decoy61_82_A_AMBER_minimized.pdb':        6.98,
'decoy61_143_A_AMBER_minimized.pdb':       8.14,
'decoy61_164_A_AMBER_minimized.pdb':       9.36,
'decoy63_1_A_AMBER_minimized.pdb':         4.80,
'decoy63_19_A_AMBER_minimized.pdb':        6.01,
'decoy65_1_A_AMBER_minimized.pdb':         6.99,
'decoy65_8_A_AMBER_minimized.pdb':         7.72,
'decoy65_11_A_AMBER_minimized.pdb':        7.86,
'decoy65_31_A_AMBER_minimized.pdb':        7.99,
'decoy65_91_A_AMBER_minimized.pdb':        8.70,
'decoy65_93_A_AMBER_minimized.pdb':        8.57,
'decoy65_52_A_AMBER_minimized.pdb':        8.62,
'decoy65_103_A_AMBER_minimized.pdb':       9.12,
'decoy65_126_A_AMBER_minimized.pdb':      10.28,
'decoy65_125_A_AMBER_minimized.pdb':      10.28,
'decoy66_1_A_AMBER_minimized.pdb':         5.55,
'decoy66_9_A_AMBER_minimized.pdb':         6.57,
'decoy66_26_A_AMBER_minimized.pdb':        7.37,
'decoy66_106_A_AMBER_minimized.pdb':       8.46,
'decoy66_112_A_AMBER_minimized.pdb':       9.81,
'decoy66_140_A_AMBER_minimized.pdb':      10.10,
'decoy67_1_A_AMBER_minimized.pdb':         5.15,
'decoy67_14_A_AMBER_minimized.pdb':        6.84,
'decoy67_147_A_AMBER_minimized.pdb':       9.96,
'decoy67_178_A_AMBER_minimized.pdb':      10.47,
'decoy68_1_A_AMBER_minimized.pdb':         5.63,
'decoy68_12_A_AMBER_minimized.pdb':        6.32,
'decoy68_15_A_AMBER_minimized.pdb':        7.91,
'decoy68_78_A_AMBER_minimized.pdb':        8.91,
'decoy68_116_A_AMBER_minimized.pdb':       9.03,
'decoy68_130_A_AMBER_minimized.pdb':       9.98,
'decoy69_1_A_AMBER_minimized.pdb':         4.81,
'decoy69_12_A_AMBER_minimized.pdb':        5.69,
'decoy69_37_A_AMBER_minimized.pdb':        6.89,
'decoy69_75_A_AMBER_minimized.pdb':        7.20,
'decoy69_131_A_AMBER_minimized.pdb':       8.64,
'decoy69_202_A_AMBER_minimized.pdb':       9.08,
'decoy69_218_A_AMBER_minimized.pdb':      10.06,
}

#print system.energy(log       = True,
#                    AB_energy = True)


#from random import random
monte_carlo(molecule           = system      ,
            #random             = random     ,
            temperature        = 100         ,
            Kb                 = 1           ,
            angle_range        = 60          ,
            nSteps             = 50000       ,
            fragment_rate      = 0.0         , #between 0  and 1
            fragment_sidechain = True        ,
            PhiPsi_rate        = 1.0         ,
            trajectory         = 'trajectory',
            pn                 = 1                       )




for decoy in decoys:
    system.load_PDB_to_system      (filename = os.path.join(folder, decoy) )   
    filein = decoy.split('/')
    filein = filein[-1]
    print  decoy, system.energy(), decoys[decoy]




'''

#-------------------------------------------------------------------------------
system = Molecule() 
system.load_PDB_to_system      (filename = os.path.join(PEPDICE_EXAMPLES , 'PDBs/alpha/1GAB/1GAB_A_AMBER_minimized.pdb'   )   )   

system.ff_type = 'Calpha_model'
system.import_Calpha_model_parameters(top       = os.path.join(PEPDICE_EXAMPLES , 'PDBs/alpha/1GAB/1GAB_A_AMBER.top')   ,   
                                      torsions  = os.path.join(PEPDICE_PARAMETER, 'amber/AMBER_rotamers.dat') )  
TRAJECTORY = 'Calpha_test'                                      
                                      

#-------------------------------------------------------------------------------
system.AB = 1
print system.energy(log       = True,
                    AB_energy = True)

'''















'''
#--------------------------------------------------------------------------------------------------------------------------------------------

files =  os.listdir('/home/fernando/programs/pepdice/Examples/teste/RO1aiu/')

rmsd  =  {}

text = open('/home/fernando/programs/pepdice/Examples/teste/RO1aiu/' + 'list.txt', 'r')
for line in text:
    line2 = line.split()
    if len(line2) == 2:
        if line2[0] == 'NAME':
            pass
        else:
            rmsd[line2[0]] = line2[1]



for _file in files:
    file2 = _file.split('.')
    
    #print file2
    
    if file2[-1] == 'pdb':
        
        system = Molecule() 
        system.load_PDB_to_system      (filename = os.path.join('/home/fernando/programs/pepdice/Examples/teste/RO1aiu/' , _file   )   )   

        system.ff_type = 'Calpha_model'
        system.import_Calpha_model_parameters(top       = None   ,   
                                              torsions  = os.path.join(PEPDICE_PARAMETER, 'amber/AMBER_rotamers.dat') )  

        print '%20s %20.10f %10s'%(_file, system.energy(log   = True,  AB_energy = True), rmsd[_file])
#--------------------------------------------------------------------------------------------------------------------------------------------
'''



'''

pdbs = [
        '/home/fernando/programs/pepdice/Examples/PDBs/alpha/1GAB/1GAB_A_AMBER.pdb']

system  = build_fragment_library_from_pdbs (
                                            molecule             = system ,
                                            frag_size            = 3      ,
                                            number_of_fragments  = 100    ,
                                            pdblist              = pdbs   ,
                                            )


import pickle

pickle.dump(system.fragments, open( "1gab_fragments.p", "wb" ) )















import pickle
system.fragments = pickle.load( open( "1gab_fragments.p", "rb" ) )

i =1
monte_carlo(molecule           = system      ,
            #random             = random      ,
            temperature        = 10000          ,
            Kb                 = 1           ,
            angle_range        = 5           ,
            nSteps             = 5000        ,
            fragment_rate      = 1.0         , #between 0  and 1
            fragment_sidechain = True        ,
            PhiPsi_rate        = 1.0         ,
            trajectory         = TRAJECTORY+str(i),
            pn                 = 1                       )

'''
