#'''
#import pandas
#from scipy.interpolate import interp2d

#filein  = '/home/farminf/Programas/pepdice/Parameters/NDRD/NDRD_TCBIG.txt'

#lines = pandas.read_csv(filein, skiprows=59, names=['res', 'dir', 'neigh', 'phi', 'psi', 'prob', 'logprob', 'cumsum'], sep='\s+')

#probabilities = {}

#for res in lines['res'].unique():
    #for neigh in ['ALL']:
        #table = lines[lines['res'] == res][lines['neigh'] == neigh][lines['dir'] == 'left']
        #probabilities[res] = interp2d(table['phi'], table['psi'], table['logprob'])

#with open('NDRD_TCBIG.txt', 'w') as output_file:
    #pickle.dump(probabilities, output_file)
    

import pickle
with open('NDRD_TCBIG.pkl', 'r') as input_file:
    probabilities = pickle.load(input_file)
    
#-------------------------------------------------------------------------------
#'''



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
        
        
        resn =  system.residues[i].name  
        
        #ala_val = lines[lines['res'] == resn][lines['neigh'] == 'ALL'][lines['dir'] == 'left']
        #print ala_val['phi'], ala_val['psi'], ala_val['logprob']
        #f = scipy.interpolate.interp2d(ala_val['phi'], ala_val['psi'], ala_val['logprob'])
        logP = probabilities[resn](phi_final_angle, psi_final_angle)
        
        print system.residues[i].name, phi_final_angle,   psi_final_angle,   logP
        torsions.append([phi_final_angle, 
                         psi_final_angle, 
                         ome_final_angle,
                         logP,
                         ])
        

    if log:
        print '''
----------------------------------------------------------
\n\n''' 
    return torsions


#eita

#-------------------------------------------------------------------------------
PEPDICE = os.environ.get('PEPDICE')
PEPDICE_EXAMPLES = os.path.join(PEPDICE, 'Examples')
PEPDICE_PARAMETER= os.path.join(PEPDICE, 'Parameters')
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
system = Molecule() 
system.load_PDB_to_system      (filename = os.path.join(PEPDICE_EXAMPLES , 'PDBs/alpha/1L2Y/1L2Y_A_AMBER_minimized.pdb'   )   )   
system.import_AMBER_parameters(top       = os.path.join(PEPDICE_EXAMPLES , 'PDBs/alpha/1L2Y/1L2Y_A_AMBER.top')   ,   
                                torsions = os.path.join(PEPDICE_PARAMETER, 'amber/AMBER_rotamers.dat') )   
#-------------------------------------------------------------------------------


pdbs = [
'1L2Y_A_AMBER_minimized.pdb'	     ,'decoy54_133_A_AMBER_minimized.pdb',
'decoy0_1_A_AMBER_minimized.pdb'	 ,'decoy54_140_A_AMBER_minimized.pdb',
'decoy10_4_A_AMBER_minimized.pdb'    ,'decoy54_163_A_AMBER_minimized.pdb'  ,
'decoy1_10_A_AMBER_minimized.pdb'    ,'decoy54_168_A_AMBER_minimized.pdb'  ,
'decoy1_11_A_AMBER_minimized.pdb'    ,'decoy54_171_A_AMBER_minimized.pdb'  ,
'decoy11_1_A_AMBER_minimized.pdb'    ,'decoy54_187_A_AMBER_minimized.pdb'  ,
'decoy1_12_A_AMBER_minimized.pdb'    ,'decoy54_191_A_AMBER_minimized.pdb'  ,
'decoy1_13_A_AMBER_minimized.pdb'    ,'decoy54_30_A_AMBER_minimized.pdb'   ,
'decoy1_14_A_AMBER_minimized.pdb'    ,'decoy54_89_A_AMBER_minimized.pdb'   ,
'decoy1_15_A_AMBER_minimized.pdb'    ,'decoy54_96_A_AMBER_minimized.pdb'   ,
'decoy1_16_A_AMBER_minimized.pdb'    ,'decoy55_19_A_AMBER_minimized.pdb'   ,
'decoy1_17_A_AMBER_minimized.pdb'    ,'decoy56_4_A_AMBER_minimized.pdb'    ,
'decoy1_18_A_AMBER_minimized.pdb'    ,'decoy57_4_A_AMBER_minimized.pdb'    ,
'decoy1_19_A_AMBER_minimized.pdb'    ,'decoy58_1_A_AMBER_minimized.pdb'    ,
'decoy1_1_A_AMBER_minimized.pdb'	   ,'decoy59_1_A_AMBER_minimized.pdb'  ,
'decoy1_20_A_AMBER_minimized.pdb'    ,'decoy60_2_A_AMBER_minimized.pdb'    ,
'decoy12_177_A_AMBER_minimized.pdb'  ,'decoy61_104_A_AMBER_minimized.pdb'  ,
'decoy1_21_A_AMBER_minimized.pdb'    ,'decoy61_115_A_AMBER_minimized.pdb'  ,
'decoy12_1_A_AMBER_minimized.pdb'    ,'decoy61_120_A_AMBER_minimized.pdb'  ,
'decoy1_22_A_AMBER_minimized.pdb'    ,'decoy61_126_A_AMBER_minimized.pdb'  ,
'decoy1_23_A_AMBER_minimized.pdb'    ,'decoy61_183_A_AMBER_minimized.pdb'  ,
'decoy1_24_A_AMBER_minimized.pdb'    ,'decoy61_4_A_AMBER_minimized.pdb'    ,
'decoy1_25_A_AMBER_minimized.pdb'    ,'decoy62_4_A_AMBER_minimized.pdb'    ,
'decoy1_26_A_AMBER_minimized.pdb'    ,'decoy6_26_A_AMBER_minimized.pdb'    ,
'decoy1_27_A_AMBER_minimized.pdb'    ,'decoy63_1_A_AMBER_minimized.pdb'    ,
'decoy1_28_A_AMBER_minimized.pdb'    ,'decoy64_119_A_AMBER_minimized.pdb'  ,
'decoy1_29_A_AMBER_minimized.pdb'    ,'decoy64_120_A_AMBER_minimized.pdb'  ,
'decoy1_2_A_AMBER_minimized.pdb'	   ,'decoy64_125_A_AMBER_minimized.pdb',
'decoy1_30_A_AMBER_minimized.pdb'    ,'decoy64_126_A_AMBER_minimized.pdb'  ,
'decoy13_10_A_AMBER_minimized.pdb'   ,'decoy64_127_A_AMBER_minimized.pdb'  ,
'decoy1_31_A_AMBER_minimized.pdb'    ,'decoy64_130_A_AMBER_minimized.pdb'  ,
'decoy1_32_A_AMBER_minimized.pdb'    ,'decoy64_131_A_AMBER_minimized.pdb'  ,
'decoy1_33_A_AMBER_minimized.pdb'    ,'decoy64_135_A_AMBER_minimized.pdb'  ,
'decoy1_34_A_AMBER_minimized.pdb'    ,'decoy64_136_A_AMBER_minimized.pdb'  ,
'decoy1_35_A_AMBER_minimized.pdb'    ,'decoy64_140_A_AMBER_minimized.pdb'  ,
'decoy1_36_A_AMBER_minimized.pdb'    ,'decoy64_142_A_AMBER_minimized.pdb'  ,
'decoy1_3_A_AMBER_minimized.pdb'	   ,'decoy64_146_A_AMBER_minimized.pdb',
'decoy14_2_A_AMBER_minimized.pdb'    ,'decoy64_147_A_AMBER_minimized.pdb'  ,
'decoy1_4_A_AMBER_minimized.pdb'	   ,'decoy64_153_A_AMBER_minimized.pdb',
'decoy15_248_A_AMBER_minimized.pdb'  ,'decoy64_154_A_AMBER_minimized.pdb'  ,
'decoy15_249_A_AMBER_minimized.pdb'  ,'decoy64_158_A_AMBER_minimized.pdb'  ,
'decoy15_250_A_AMBER_minimized.pdb'  ,'decoy64_159_A_AMBER_minimized.pdb'  ,
'decoy15_251_A_AMBER_minimized.pdb'  ,'decoy64_165_A_AMBER_minimized.pdb'  ,
'decoy15_40_A_AMBER_minimized.pdb'   ,'decoy64_169_A_AMBER_minimized.pdb'  ,
'decoy1_5_A_AMBER_minimized.pdb'	   ,'decoy64_170_A_AMBER_minimized.pdb',
'decoy16_1_A_AMBER_minimized.pdb'    ,'decoy64_172_A_AMBER_minimized.pdb'  ,
'decoy1_6_A_AMBER_minimized.pdb'	   ,'decoy64_177_A_AMBER_minimized.pdb',
'decoy17_224_A_AMBER_minimized.pdb'  ,'decoy64_178_A_AMBER_minimized.pdb'  ,
'decoy17_225_A_AMBER_minimized.pdb'  ,'decoy64_179_A_AMBER_minimized.pdb'  ,
'decoy17_226_A_AMBER_minimized.pdb'  ,'decoy64_180_A_AMBER_minimized.pdb'  ,
'decoy17_227_A_AMBER_minimized.pdb'  ,'decoy64_182_A_AMBER_minimized.pdb'  ,
'decoy17_228_A_AMBER_minimized.pdb'  ,'decoy64_183_A_AMBER_minimized.pdb'  ,
'decoy17_229_A_AMBER_minimized.pdb'  ,'decoy64_184_A_AMBER_minimized.pdb'  ,
'decoy17_48_A_AMBER_minimized.pdb'   ,'decoy64_185_A_AMBER_minimized.pdb'  ,
'decoy1_7_A_AMBER_minimized.pdb'	   ,'decoy64_188_A_AMBER_minimized.pdb',
'decoy18_8_A_AMBER_minimized.pdb'    ,'decoy64_189_A_AMBER_minimized.pdb'  ,
'decoy1_8_A_AMBER_minimized.pdb'	   ,'decoy64_191_A_AMBER_minimized.pdb',
'decoy19_30_A_AMBER_minimized.pdb'   ,'decoy64_192_A_AMBER_minimized.pdb'  ,
'decoy1_9_A_AMBER_minimized.pdb'	   ,'decoy64_193_A_AMBER_minimized.pdb',
'decoy20_45_A_AMBER_minimized.pdb'   ,'decoy64_194_A_AMBER_minimized.pdb'  ,
'decoy2_11_A_AMBER_minimized.pdb'    ,'decoy64_195_A_AMBER_minimized.pdb'  ,
'decoy21_32_A_AMBER_minimized.pdb'   ,'decoy64_196_A_AMBER_minimized.pdb'  ,
'decoy22_31_A_AMBER_minimized.pdb'   ,'decoy64_197_A_AMBER_minimized.pdb'  ,
'decoy2_249_A_AMBER_minimized.pdb'   ,'decoy64_198_A_AMBER_minimized.pdb'  ,
'decoy2_251_A_AMBER_minimized.pdb'   ,'decoy64_199_A_AMBER_minimized.pdb'  ,
'decoy23_1_A_AMBER_minimized.pdb'    ,'decoy64_1_A_AMBER_minimized.pdb'    ,
'decoy24_15_A_AMBER_minimized.pdb'   ,'decoy64_200_A_AMBER_minimized.pdb'  ,
'decoy25_35_A_AMBER_minimized.pdb'   ,'decoy64_201_A_AMBER_minimized.pdb'  ,
'decoy26_104_A_AMBER_minimized.pdb'  ,'decoy64_202_A_AMBER_minimized.pdb'  ,
'decoy26_152_A_AMBER_minimized.pdb'  ,'decoy64_203_A_AMBER_minimized.pdb'  ,
'decoy26_2_A_AMBER_minimized.pdb'    ,'decoy64_204_A_AMBER_minimized.pdb'  ,
'decoy27_164_A_AMBER_minimized.pdb'  ,'decoy65_20_A_AMBER_minimized.pdb'   ,
'decoy27_201_A_AMBER_minimized.pdb'  ,'decoy66_14_A_AMBER_minimized.pdb'   ,
'decoy27_50_A_AMBER_minimized.pdb'   ,'decoy66_194_A_AMBER_minimized.pdb'  ,
'decoy28_160_A_AMBER_minimized.pdb'  ,'decoy66_195_A_AMBER_minimized.pdb'  ,
'decoy28_1_A_AMBER_minimized.pdb'    ,'decoy66_196_A_AMBER_minimized.pdb'  ,
'decoy29_21_A_AMBER_minimized.pdb'   ,'decoy66_197_A_AMBER_minimized.pdb'  ,
'decoy30_56_A_AMBER_minimized.pdb'   ,'decoy67_108_A_AMBER_minimized.pdb'  ,
'decoy31_8_A_AMBER_minimized.pdb'    ,'decoy67_141_A_AMBER_minimized.pdb'  ,
'decoy3_1_A_AMBER_minimized.pdb'	   ,'decoy67_145_A_AMBER_minimized.pdb',
'decoy32_21_A_AMBER_minimized.pdb'   ,'decoy67_147_A_AMBER_minimized.pdb'  ,
'decoy33_26_A_AMBER_minimized.pdb'   ,'decoy67_162_A_AMBER_minimized.pdb'  ,
'decoy34_1_A_AMBER_minimized.pdb'    ,'decoy67_165_A_AMBER_minimized.pdb'  ,
'decoy35_42_A_AMBER_minimized.pdb'   ,'decoy67_177_A_AMBER_minimized.pdb'  ,
'decoy35_92_A_AMBER_minimized.pdb'   ,'decoy67_178_A_AMBER_minimized.pdb'  ,
'decoy36_15_A_AMBER_minimized.pdb'   ,'decoy67_179_A_AMBER_minimized.pdb'  ,
'decoy36_233_A_AMBER_minimized.pdb'  ,'decoy67_180_A_AMBER_minimized.pdb'  ,
'decoy36_87_A_AMBER_minimized.pdb'   ,'decoy67_181_A_AMBER_minimized.pdb'  ,
'decoy37_43_A_AMBER_minimized.pdb'   ,'decoy67_198_A_AMBER_minimized.pdb'  ,
'decoy38_28_A_AMBER_minimized.pdb'   ,'decoy67_199_A_AMBER_minimized.pdb'  ,
'decoy39_23_A_AMBER_minimized.pdb'   ,'decoy67_200_A_AMBER_minimized.pdb'  ,
'decoy40_152_A_AMBER_minimized.pdb'  ,'decoy67_26_A_AMBER_minimized.pdb'   ,
'decoy40_153_A_AMBER_minimized.pdb'  ,'decoy68_1_A_AMBER_minimized.pdb'    ,
'decoy40_1_A_AMBER_minimized.pdb'    ,'decoy69_100_A_AMBER_minimized.pdb'  ,
'decoy4_125_A_AMBER_minimized.pdb'   ,'decoy69_102_A_AMBER_minimized.pdb'  ,
'decoy41_8_A_AMBER_minimized.pdb'    ,'decoy69_103_A_AMBER_minimized.pdb'  ,
'decoy4_201_A_AMBER_minimized.pdb'   ,'decoy69_105_A_AMBER_minimized.pdb'  ,
'decoy4_224_A_AMBER_minimized.pdb'   ,'decoy69_106_A_AMBER_minimized.pdb'  ,
'decoy42_51_A_AMBER_minimized.pdb'   ,'decoy69_108_A_AMBER_minimized.pdb'  ,
'decoy4_262_A_AMBER_minimized.pdb'   ,'decoy69_111_A_AMBER_minimized.pdb'  ,
'decoy4_270_A_AMBER_minimized.pdb'   ,'decoy69_112_A_AMBER_minimized.pdb'  ,
'decoy43_3_A_AMBER_minimized.pdb'    ,'decoy69_122_A_AMBER_minimized.pdb'  ,
'decoy44_50_A_AMBER_minimized.pdb'   ,'decoy69_151_A_AMBER_minimized.pdb'  ,
'decoy45_101_A_AMBER_minimized.pdb'  ,'decoy69_154_A_AMBER_minimized.pdb'  ,
'decoy45_104_A_AMBER_minimized.pdb'  ,'decoy69_155_A_AMBER_minimized.pdb'  ,
'decoy45_160_A_AMBER_minimized.pdb'  ,'decoy69_159_A_AMBER_minimized.pdb'  ,
'decoy45_1_A_AMBER_minimized.pdb'    ,'decoy69_160_A_AMBER_minimized.pdb'  ,
'decoy46_1_A_AMBER_minimized.pdb'    ,'decoy69_193_A_AMBER_minimized.pdb'  ,
'decoy47_1_A_AMBER_minimized.pdb'    ,'decoy69_194_A_AMBER_minimized.pdb'  ,
'decoy48_42_A_AMBER_minimized.pdb'   ,'decoy69_196_A_AMBER_minimized.pdb'  ,
'decoy49_2_A_AMBER_minimized.pdb'    ,'decoy69_243_A_AMBER_minimized.pdb'  ,
'decoy49_58_A_AMBER_minimized.pdb'   ,'decoy69_245_A_AMBER_minimized.pdb'  ,
'decoy49_59_A_AMBER_minimized.pdb'   ,'decoy69_246_A_AMBER_minimized.pdb'  ,
'decoy49_64_A_AMBER_minimized.pdb'   ,'decoy69_38_A_AMBER_minimized.pdb'   ,
'decoy49_69_A_AMBER_minimized.pdb'   ,'decoy69_95_A_AMBER_minimized.pdb'   ,
'decoy49_72_A_AMBER_minimized.pdb'   ,'decoy69_96_A_AMBER_minimized.pdb'   ,
'decoy49_74_A_AMBER_minimized.pdb'   ,'decoy69_98_A_AMBER_minimized.pdb'   ,
'decoy49_87_A_AMBER_minimized.pdb'   ,'decoy70_13_A_AMBER_minimized.pdb'   ,
'decoy49_90_A_AMBER_minimized.pdb'   ,'decoy71_27_A_AMBER_minimized.pdb'   ,
'decoy4_9_A_AMBER_minimized.pdb'	   ,'decoy7_1_A_AMBER_minimized.pdb'   ,
'decoy50_176_A_AMBER_minimized.pdb'  ,'decoy72_34_A_AMBER_minimized.pdb'   ,
'decoy50_177_A_AMBER_minimized.pdb'  ,'decoy73_16_A_AMBER_minimized.pdb'   ,
'decoy50_3_A_AMBER_minimized.pdb'    ,'decoy73_207_A_AMBER_minimized.pdb'  ,
'decoy51_56_A_AMBER_minimized.pdb'   ,'decoy73_208_A_AMBER_minimized.pdb'  ,
'decoy52_1_A_AMBER_minimized.pdb'    ,'decoy73_210_A_AMBER_minimized.pdb'  ,
'decoy5_2_A_AMBER_minimized.pdb'	   ,'decoy73_211_A_AMBER_minimized.pdb',
'decoy53_161_A_AMBER_minimized.pdb'  ,'decoy73_212_A_AMBER_minimized.pdb'  ,
'decoy53_163_A_AMBER_minimized.pdb'  ,'decoy73_213_A_AMBER_minimized.pdb'  ,
'decoy53_168_A_AMBER_minimized.pdb'  ,'decoy73_214_A_AMBER_minimized.pdb'  ,
'decoy53_169_A_AMBER_minimized.pdb'  ,'decoy73_215_A_AMBER_minimized.pdb'  ,
'decoy53_171_A_AMBER_minimized.pdb'  ,'decoy74_14_A_AMBER_minimized.pdb'   ,
'decoy53_177_A_AMBER_minimized.pdb'  ,'decoy75_3_A_AMBER_minimized.pdb'    ,
'decoy53_181_A_AMBER_minimized.pdb'  ,'decoy76_2_A_AMBER_minimized.pdb'    ,
'decoy53_182_A_AMBER_minimized.pdb'  ,'decoy77_18_A_AMBER_minimized.pdb'   ,
'decoy53_184_A_AMBER_minimized.pdb'  ,'decoy78_5_A_AMBER_minimized.pdb'    ,
'decoy53_185_A_AMBER_minimized.pdb'  ,'decoy79_136_A_AMBER_minimized.pdb'  ,
'decoy53_186_A_AMBER_minimized.pdb'  ,'decoy79_2_A_AMBER_minimized.pdb'    ,
'decoy53_189_A_AMBER_minimized.pdb'  ,'decoy80_7_A_AMBER_minimized.pdb'    ,
'decoy53_191_A_AMBER_minimized.pdb'  ,'decoy81_89_A_AMBER_minimized.pdb'   ,
'decoy53_193_A_AMBER_minimized.pdb'  ,'decoy8_1_A_AMBER_minimized.pdb'     ,
'decoy53_197_A_AMBER_minimized.pdb'  ,'decoy82_191_A_AMBER_minimized.pdb'  ,
'decoy53_198_A_AMBER_minimized.pdb'  ,'decoy82_27_A_AMBER_minimized.pdb'   ,
'decoy53_199_A_AMBER_minimized.pdb'  ,'decoy83_54_A_AMBER_minimized.pdb'   ,
'decoy53_203_A_AMBER_minimized.pdb'  ,'decoy84_14_A_AMBER_minimized.pdb'   ,
'decoy53_211_A_AMBER_minimized.pdb'  ,'decoy85_52_A_AMBER_minimized.pdb'   ,
'decoy53_5_A_AMBER_minimized.pdb'    ,'decoy86_2_A_AMBER_minimized.pdb'    ,
'decoy54_103_A_AMBER_minimized.pdb'  ,'decoy87_26_A_AMBER_minimized.pdb'   ,
'decoy54_106_A_AMBER_minimized.pdb'  ,'decoy9_123_A_AMBER_minimized.pdb'   ,
'decoy54_112_A_AMBER_minimized.pdb'  ,'decoy9_143_A_AMBER_minimized.pdb'   ,
'decoy54_115_A_AMBER_minimized.pdb'  ,'decoy9_18_A_AMBER_minimized.pdb'    ,
'decoy54_118_A_AMBER_minimized.pdb']

rmsd = {
'1L2Y_A_AMBER_minimized.pdb'           : 0.00,
'decoy0_1_A_AMBER_minimized.pdb'       : 0.05,
'decoy1_2_A_AMBER_minimized.pdb'       : 1.08,
'decoy1_8_A_AMBER_minimized.pdb'       : 0.81,
'decoy1_9_A_AMBER_minimized.pdb'       : 1.02,
'decoy1_11_A_AMBER_minimized.pdb'      : 0.81,
'decoy1_12_A_AMBER_minimized.pdb'      : 0.86,
'decoy1_13_A_AMBER_minimized.pdb'      : 0.97,
'decoy1_22_A_AMBER_minimized.pdb'      : 0.94,
'decoy1_24_A_AMBER_minimized.pdb'      : 0.82,
'decoy1_26_A_AMBER_minimized.pdb'      : 0.77,
'decoy1_30_A_AMBER_minimized.pdb'      : 0.91,
'decoy1_34_A_AMBER_minimized.pdb'      : 0.91,
'decoy1_3_A_AMBER_minimized.pdb'       : 1.05,
'decoy1_4_A_AMBER_minimized.pdb'       : 0.84,
'decoy1_5_A_AMBER_minimized.pdb'       : 0.84,
'decoy1_6_A_AMBER_minimized.pdb'       : 0.49,
'decoy1_10_A_AMBER_minimized.pdb'      : 0.56,
'decoy1_16_A_AMBER_minimized.pdb'      : 0.93,
'decoy1_17_A_AMBER_minimized.pdb'      : 0.98,
'decoy1_19_A_AMBER_minimized.pdb'      : 0.93,
'decoy1_20_A_AMBER_minimized.pdb'      : 0.62,
'decoy1_21_A_AMBER_minimized.pdb'      : 0.74,
'decoy1_23_A_AMBER_minimized.pdb'      : 0.54,
'decoy1_25_A_AMBER_minimized.pdb'      : 0.81,
'decoy1_27_A_AMBER_minimized.pdb'      : 0.72,
'decoy1_29_A_AMBER_minimized.pdb'      : 1.02,
'decoy1_31_A_AMBER_minimized.pdb'      : 0.90,
'decoy1_32_A_AMBER_minimized.pdb'      : 0.97,
'decoy1_35_A_AMBER_minimized.pdb'      : 0.74,
'decoy1_7_A_AMBER_minimized.pdb'       : 0.72,
'decoy1_14_A_AMBER_minimized.pdb'      : 0.58,
'decoy1_15_A_AMBER_minimized.pdb'      : 0.76,
'decoy1_18_A_AMBER_minimized.pdb'      : 0.81,
'decoy1_28_A_AMBER_minimized.pdb'      : 0.58,
'decoy1_1_A_AMBER_minimized.pdb'       : 1.07,
'decoy1_33_A_AMBER_minimized.pdb'      : 0.97,
'decoy1_36_A_AMBER_minimized.pdb'      : 0.90,
'decoy2_11_A_AMBER_minimized.pdb'      : 1.95,
'decoy2_249_A_AMBER_minimized.pdb'     : 6.28,
'decoy2_251_A_AMBER_minimized.pdb'     : 6.27,
'decoy3_1_A_AMBER_minimized.pdb'       : 1.90,
'decoy4_9_A_AMBER_minimized.pdb'       : 1.87,
'decoy4_125_A_AMBER_minimized.pdb'     : 3.87,
'decoy4_201_A_AMBER_minimized.pdb'     : 5.12,
'decoy4_224_A_AMBER_minimized.pdb'     : 5.78,
'decoy4_262_A_AMBER_minimized.pdb'     : 6.08,
'decoy4_270_A_AMBER_minimized.pdb'     : 7.45,
'decoy5_2_A_AMBER_minimized.pdb'       : 1.60,
'decoy6_26_A_AMBER_minimized.pdb'      : 1.91,
'decoy7_1_A_AMBER_minimized.pdb'       : 1.98,
'decoy8_1_A_AMBER_minimized.pdb'       : 1.91,
'decoy9_18_A_AMBER_minimized.pdb'      : 1.89,
'decoy9_123_A_AMBER_minimized.pdb'     : 3.65,
'decoy9_143_A_AMBER_minimized.pdb'     : 3.80,
'decoy10_4_A_AMBER_minimized.pdb'      : 1.74,
'decoy11_1_A_AMBER_minimized.pdb'      : 1.91,
'decoy12_1_A_AMBER_minimized.pdb'      : 1.85,
'decoy12_177_A_AMBER_minimized.pdb'    : 7.24,
'decoy13_10_A_AMBER_minimized.pdb'     : 1.53,
'decoy14_2_A_AMBER_minimized.pdb'      : 1.54,
'decoy15_40_A_AMBER_minimized.pdb'     : 1.94,
'decoy15_248_A_AMBER_minimized.pdb'    : 7.62,
'decoy15_249_A_AMBER_minimized.pdb'    : 7.99,
'decoy15_250_A_AMBER_minimized.pdb'    : 7.23,
'decoy15_251_A_AMBER_minimized.pdb'    : 8.00,
'decoy16_1_A_AMBER_minimized.pdb'      : 1.99,
'decoy17_48_A_AMBER_minimized.pdb'     : 1.69,
'decoy17_228_A_AMBER_minimized.pdb'    : 7.20,
'decoy17_229_A_AMBER_minimized.pdb'    : 7.03,
'decoy17_224_A_AMBER_minimized.pdb'    : 7.07,
'decoy17_225_A_AMBER_minimized.pdb'    : 7.40,
'decoy17_226_A_AMBER_minimized.pdb'    : 6.90,
'decoy17_227_A_AMBER_minimized.pdb'    : 7.12,
'decoy18_8_A_AMBER_minimized.pdb'      : 1.75,
'decoy19_30_A_AMBER_minimized.pdb'     : 1.63,
'decoy20_45_A_AMBER_minimized.pdb'     : 1.80,
'decoy21_32_A_AMBER_minimized.pdb'     : 1.90,
'decoy22_31_A_AMBER_minimized.pdb'     : 1.86,
'decoy23_1_A_AMBER_minimized.pdb'      : 1.92,
'decoy24_15_A_AMBER_minimized.pdb'     : 1.95,
'decoy25_35_A_AMBER_minimized.pdb'     : 1.83,
'decoy26_2_A_AMBER_minimized.pdb'      : 2.93,
'decoy26_104_A_AMBER_minimized.pdb'    : 4.43,
'decoy26_152_A_AMBER_minimized.pdb'    : 5.17,
'decoy27_50_A_AMBER_minimized.pdb'     : 1.87,
'decoy27_164_A_AMBER_minimized.pdb'    : 4.70,
'decoy27_201_A_AMBER_minimized.pdb'    : 5.15,
'decoy28_1_A_AMBER_minimized.pdb'      : 1.88,
'decoy28_160_A_AMBER_minimized.pdb'    : 6.30,
'decoy29_21_A_AMBER_minimized.pdb'     : 1.90,
'decoy30_56_A_AMBER_minimized.pdb'     : 2.60,
'decoy31_8_A_AMBER_minimized.pdb'      : 1.96,
'decoy32_21_A_AMBER_minimized.pdb'     : 1.73,
'decoy33_26_A_AMBER_minimized.pdb'     : 2.02,
'decoy34_1_A_AMBER_minimized.pdb'      : 1.86,
'decoy35_42_A_AMBER_minimized.pdb'     : 2.65,
'decoy35_92_A_AMBER_minimized.pdb'     : 5.58,
'decoy36_15_A_AMBER_minimized.pdb'     : 1.95,
'decoy36_87_A_AMBER_minimized.pdb'     : 4.00,
'decoy36_233_A_AMBER_minimized.pdb'    : 6.40,
'decoy37_43_A_AMBER_minimized.pdb'     : 2.71,
'decoy38_28_A_AMBER_minimized.pdb'     : 2.66,
'decoy39_23_A_AMBER_minimized.pdb'     : 1.94,
'decoy40_1_A_AMBER_minimized.pdb'      : 2.95,
'decoy40_152_A_AMBER_minimized.pdb'    : 5.40,
'decoy40_153_A_AMBER_minimized.pdb'    : 5.48,
'decoy41_8_A_AMBER_minimized.pdb'      : 2.58,
'decoy42_51_A_AMBER_minimized.pdb'     : 3.95,
'decoy43_3_A_AMBER_minimized.pdb'      : 1.75,
'decoy44_50_A_AMBER_minimized.pdb'     : 2.73,
'decoy45_1_A_AMBER_minimized.pdb'      : 2.85,
'decoy45_101_A_AMBER_minimized.pdb'    : 4.35,
'decoy45_104_A_AMBER_minimized.pdb'    : 4.51,
'decoy45_160_A_AMBER_minimized.pdb'    : 6.62,
'decoy46_1_A_AMBER_minimized.pdb'      : 3.02,
'decoy47_1_A_AMBER_minimized.pdb'      : 3.02,
'decoy48_42_A_AMBER_minimized.pdb'     : 2.80,
'decoy49_2_A_AMBER_minimized.pdb'      : 2.78,
'decoy49_58_A_AMBER_minimized.pdb'     : 3.41,
'decoy49_59_A_AMBER_minimized.pdb'     : 3.38,
'decoy49_64_A_AMBER_minimized.pdb'     : 3.35,
'decoy49_69_A_AMBER_minimized.pdb'     : 3.59,
'decoy49_72_A_AMBER_minimized.pdb'     : 3.87,
'decoy49_74_A_AMBER_minimized.pdb'     : 3.73,
'decoy49_87_A_AMBER_minimized.pdb'     : 3.44,
'decoy49_90_A_AMBER_minimized.pdb'     : 3.94,
'decoy50_3_A_AMBER_minimized.pdb'      : 2.86,
'decoy50_176_A_AMBER_minimized.pdb'    : 7.05,
'decoy50_177_A_AMBER_minimized.pdb'    : 7.19,
'decoy51_56_A_AMBER_minimized.pdb'     : 3.67,
'decoy52_1_A_AMBER_minimized.pdb'      : 2.95,
'decoy53_5_A_AMBER_minimized.pdb'      : 2.99,
'decoy53_177_A_AMBER_minimized.pdb'    : 6.50,
'decoy53_181_A_AMBER_minimized.pdb'    : 6.10,
'decoy53_186_A_AMBER_minimized.pdb'    : 6.67,
'decoy53_191_A_AMBER_minimized.pdb'    : 6.63,
'decoy53_193_A_AMBER_minimized.pdb'    : 6.15,
'decoy53_197_A_AMBER_minimized.pdb'    : 6.73,
'decoy53_203_A_AMBER_minimized.pdb'    : 6.60,
'decoy53_161_A_AMBER_minimized.pdb'    : 6.35,
'decoy53_163_A_AMBER_minimized.pdb'    : 6.33,
'decoy53_184_A_AMBER_minimized.pdb'    : 6.09,
'decoy53_185_A_AMBER_minimized.pdb'    : 6.10,
'decoy53_198_A_AMBER_minimized.pdb'    : 6.19,
'decoy53_199_A_AMBER_minimized.pdb'    : 6.41,
'decoy53_168_A_AMBER_minimized.pdb'    : 6.23,
'decoy53_169_A_AMBER_minimized.pdb'    : 6.53,
'decoy53_171_A_AMBER_minimized.pdb'    : 6.11,
'decoy53_182_A_AMBER_minimized.pdb'    : 6.64,
'decoy53_189_A_AMBER_minimized.pdb'    : 6.25,
'decoy53_211_A_AMBER_minimized.pdb'    : 7.00,
'decoy54_30_A_AMBER_minimized.pdb'     : 1.69,
'decoy54_89_A_AMBER_minimized.pdb'     : 3.98,
'decoy54_96_A_AMBER_minimized.pdb'     : 3.56,
'decoy54_103_A_AMBER_minimized.pdb'    : 3.51,
'decoy54_106_A_AMBER_minimized.pdb'    : 3.77,
'decoy54_112_A_AMBER_minimized.pdb'    : 3.78,
'decoy54_115_A_AMBER_minimized.pdb'    : 3.71,
'decoy54_118_A_AMBER_minimized.pdb'    : 3.55,
'decoy54_133_A_AMBER_minimized.pdb'    : 3.53,
'decoy54_140_A_AMBER_minimized.pdb'    : 4.51,
'decoy54_163_A_AMBER_minimized.pdb'    : 4.07,
'decoy54_168_A_AMBER_minimized.pdb'    : 4.53,
'decoy54_171_A_AMBER_minimized.pdb'    : 4.28,
'decoy54_187_A_AMBER_minimized.pdb'    : 4.05,
'decoy54_191_A_AMBER_minimized.pdb'    : 5.19,
'decoy55_19_A_AMBER_minimized.pdb'     : 2.81,
'decoy56_4_A_AMBER_minimized.pdb'      : 1.97,
'decoy57_4_A_AMBER_minimized.pdb'      : 2.85,
'decoy58_1_A_AMBER_minimized.pdb'      : 2.76,
'decoy59_1_A_AMBER_minimized.pdb'      : 3.02,
'decoy60_2_A_AMBER_minimized.pdb'      : 3.01,
'decoy61_4_A_AMBER_minimized.pdb'      : 1.89,
'decoy61_104_A_AMBER_minimized.pdb'    : 3.96,
'decoy61_115_A_AMBER_minimized.pdb'    : 4.46,
'decoy61_120_A_AMBER_minimized.pdb'    : 4.20,
'decoy61_126_A_AMBER_minimized.pdb'    : 4.47,
'decoy61_183_A_AMBER_minimized.pdb'    : 5.08,
'decoy62_4_A_AMBER_minimized.pdb'      : 3.04,
'decoy63_1_A_AMBER_minimized.pdb'      : 2.91,
'decoy64_1_A_AMBER_minimized.pdb'      : 2.71,
'decoy64_119_A_AMBER_minimized.pdb'    : 4.78,
'decoy64_120_A_AMBER_minimized.pdb'    : 4.87,
'decoy64_125_A_AMBER_minimized.pdb'    : 4.73,
'decoy64_126_A_AMBER_minimized.pdb'    : 4.99,
'decoy64_127_A_AMBER_minimized.pdb'    : 4.76,
'decoy64_130_A_AMBER_minimized.pdb'    : 4.89,
'decoy64_131_A_AMBER_minimized.pdb'    : 4.77,
'decoy64_135_A_AMBER_minimized.pdb'    : 4.89,
'decoy64_136_A_AMBER_minimized.pdb'    : 4.41,
'decoy64_140_A_AMBER_minimized.pdb'    : 4.63,
'decoy64_142_A_AMBER_minimized.pdb'    : 4.96,
'decoy64_146_A_AMBER_minimized.pdb'    : 4.83,
'decoy64_147_A_AMBER_minimized.pdb'    : 4.89,
'decoy64_153_A_AMBER_minimized.pdb'    : 5.31,
'decoy64_154_A_AMBER_minimized.pdb'    : 5.21,
'decoy64_165_A_AMBER_minimized.pdb'    : 5.48,
'decoy64_169_A_AMBER_minimized.pdb'    : 5.16,
'decoy64_170_A_AMBER_minimized.pdb'    : 5.29,
'decoy64_177_A_AMBER_minimized.pdb'    : 5.21,
'decoy64_184_A_AMBER_minimized.pdb'    : 5.12,
'decoy64_188_A_AMBER_minimized.pdb'    : 5.11,
'decoy64_158_A_AMBER_minimized.pdb'    : 5.86,
'decoy64_159_A_AMBER_minimized.pdb'    : 5.12,
'decoy64_172_A_AMBER_minimized.pdb'    : 5.62,
'decoy64_178_A_AMBER_minimized.pdb'    : 5.02,
'decoy64_179_A_AMBER_minimized.pdb'    : 5.60,
'decoy64_180_A_AMBER_minimized.pdb'    : 5.74,
'decoy64_182_A_AMBER_minimized.pdb'    : 5.50,
'decoy64_183_A_AMBER_minimized.pdb'    : 5.50,
'decoy64_185_A_AMBER_minimized.pdb'    : 5.00,
'decoy64_189_A_AMBER_minimized.pdb'    : 5.47,
'decoy64_191_A_AMBER_minimized.pdb'    : 5.10,
'decoy64_192_A_AMBER_minimized.pdb'    : 5.32,
'decoy64_193_A_AMBER_minimized.pdb'    : 5.34,
'decoy64_194_A_AMBER_minimized.pdb'    : 5.29,
'decoy64_195_A_AMBER_minimized.pdb'    : 5.25,
'decoy64_196_A_AMBER_minimized.pdb'    : 5.70,
'decoy64_197_A_AMBER_minimized.pdb'    : 5.29,
'decoy64_198_A_AMBER_minimized.pdb'    : 5.54,
'decoy64_199_A_AMBER_minimized.pdb'    : 5.14,
'decoy64_200_A_AMBER_minimized.pdb'    : 5.17,
'decoy64_201_A_AMBER_minimized.pdb'    : 6.23,
'decoy64_202_A_AMBER_minimized.pdb'    : 6.29,
'decoy64_203_A_AMBER_minimized.pdb'    : 6.33,
'decoy64_204_A_AMBER_minimized.pdb'    : 6.33,
'decoy65_20_A_AMBER_minimized.pdb'     : 3.06,
'decoy66_14_A_AMBER_minimized.pdb'     : 2.96,
'decoy66_194_A_AMBER_minimized.pdb'    : 7.47,
'decoy66_195_A_AMBER_minimized.pdb'    : 7.03,
'decoy66_196_A_AMBER_minimized.pdb'    : 7.35,
'decoy66_197_A_AMBER_minimized.pdb'    : 7.41,
'decoy67_26_A_AMBER_minimized.pdb'     : 1.91,
'decoy67_108_A_AMBER_minimized.pdb'    : 3.75,
'decoy67_141_A_AMBER_minimized.pdb'    : 4.67,
'decoy67_145_A_AMBER_minimized.pdb'    : 4.47,
'decoy67_147_A_AMBER_minimized.pdb'    : 4.53,
'decoy67_162_A_AMBER_minimized.pdb'    : 4.33,
'decoy67_165_A_AMBER_minimized.pdb'    : 4.39,
'decoy67_177_A_AMBER_minimized.pdb'    : 4.75,
'decoy67_178_A_AMBER_minimized.pdb'    : 4.90,
'decoy67_179_A_AMBER_minimized.pdb'    : 4.76,
'decoy67_180_A_AMBER_minimized.pdb'    : 4.91,
'decoy67_181_A_AMBER_minimized.pdb'    : 4.54,
'decoy67_198_A_AMBER_minimized.pdb'    : 6.70,
'decoy67_199_A_AMBER_minimized.pdb'    : 7.12,
'decoy67_200_A_AMBER_minimized.pdb'    : 7.29,
'decoy68_1_A_AMBER_minimized.pdb'      : 2.90,
'decoy69_38_A_AMBER_minimized.pdb'     : 1.75,
'decoy69_95_A_AMBER_minimized.pdb'     : 3.92,
'decoy69_96_A_AMBER_minimized.pdb'     : 3.87,
'decoy69_98_A_AMBER_minimized.pdb'     : 3.82,
'decoy69_100_A_AMBER_minimized.pdb'    : 3.89,
'decoy69_102_A_AMBER_minimized.pdb'    : 3.86,
'decoy69_103_A_AMBER_minimized.pdb'    : 3.79,
'decoy69_105_A_AMBER_minimized.pdb'    : 3.92,
'decoy69_106_A_AMBER_minimized.pdb'    : 3.91,
'decoy69_108_A_AMBER_minimized.pdb'    : 3.86,
'decoy69_111_A_AMBER_minimized.pdb'    : 3.92,
'decoy69_112_A_AMBER_minimized.pdb'    : 3.87,
'decoy69_122_A_AMBER_minimized.pdb'    : 3.97,
'decoy69_151_A_AMBER_minimized.pdb'    : 4.18,
'decoy69_154_A_AMBER_minimized.pdb'    : 3.93,
'decoy69_155_A_AMBER_minimized.pdb'    : 3.91,
'decoy69_159_A_AMBER_minimized.pdb'    : 3.98,
'decoy69_160_A_AMBER_minimized.pdb'    : 3.93,
'decoy69_196_A_AMBER_minimized.pdb'    : 5.52,
'decoy69_193_A_AMBER_minimized.pdb'    : 5.70,
'decoy69_194_A_AMBER_minimized.pdb'    : 5.44,
'decoy69_243_A_AMBER_minimized.pdb'    : 6.21,
'decoy69_245_A_AMBER_minimized.pdb'    : 6.03,
'decoy69_246_A_AMBER_minimized.pdb'    : 6.26,
'decoy70_13_A_AMBER_minimized.pdb'     : 2.88,
'decoy71_27_A_AMBER_minimized.pdb'     : 2.89,
'decoy72_34_A_AMBER_minimized.pdb'     : 2.94,
'decoy73_16_A_AMBER_minimized.pdb'     : 3.03,
'decoy73_207_A_AMBER_minimized.pdb'    : 6.10,
'decoy73_208_A_AMBER_minimized.pdb'    : 6.24,
'decoy73_210_A_AMBER_minimized.pdb'    : 6.00,
'decoy73_211_A_AMBER_minimized.pdb'    : 6.09,
'decoy73_214_A_AMBER_minimized.pdb'    : 5.97,
'decoy73_212_A_AMBER_minimized.pdb'    : 5.99,
'decoy73_213_A_AMBER_minimized.pdb'    : 5.94,
'decoy73_215_A_AMBER_minimized.pdb'    : 7.07,
'decoy74_14_A_AMBER_minimized.pdb'     : 2.06,
'decoy75_3_A_AMBER_minimized.pdb'      : 3.02,
'decoy76_2_A_AMBER_minimized.pdb'      : 2.92,
'decoy77_18_A_AMBER_minimized.pdb'     : 2.94,
'decoy78_5_A_AMBER_minimized.pdb'      : 2.94,
'decoy79_2_A_AMBER_minimized.pdb'      : 2.51,
'decoy79_136_A_AMBER_minimized.pdb'    : 4.43,
'decoy80_7_A_AMBER_minimized.pdb'      : 3.79,
'decoy81_89_A_AMBER_minimized.pdb'     : 3.58,
'decoy82_27_A_AMBER_minimized.pdb'     : 2.93,
'decoy82_191_A_AMBER_minimized.pdb'    : 6.94,
'decoy83_54_A_AMBER_minimized.pdb'     : 2.98,
'decoy84_14_A_AMBER_minimized.pdb'     : 2.70,
'decoy85_52_A_AMBER_minimized.pdb'     : 3.51,
'decoy86_2_A_AMBER_minimized.pdb'      : 2.67,
'decoy87_26_A_AMBER_minimized.pdb'     : 3.92,
}








text = ''

for pdb in pdbs:
    system.load_PDB_to_system      (filename = os.path.join(PEPDICE_EXAMPLES , 'PDBs/alpha/1L2Y/', pdb   )   )   
    energy = system.energy(log       = True,
                           AB_energy = True)
    
    text +=  pdb + '  ' + str(energy) + '  ' + str(rmsd[pdb]) + '\n'
                 

logfile =  open('logfile', 'w')
logfile.write(text)
#torsions = compute_torsions (system = system, log =False)

#pprint(torsions)

#system.AB = 1
#print system.energy(log       = True,
                    #AB_energy = True)
