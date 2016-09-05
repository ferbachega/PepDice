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

from MonteCarlo  import monte_carlo #,  monte_carlo_dic, MC_replica_exchange, run_MC_replica_exchange #

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




#-----------------------------------------------------------------------------------------------------------------------------------#
system = Molecule()                                                                                                                 #
system.load_PDB_to_system      (filename = os.path.join(PEPDICE_EXAMPLES , 'data/alpha/1GAB/1gab_ff03ua_AMBER_folded.pdb'     )   ) #
system.import_AMBER_parameters(top       = os.path.join(PEPDICE_EXAMPLES , 'data/alpha/1GAB/1gab_ff03ua_AMBER_extended.prmtop')   , #
                                torsions = os.path.join(PEPDICE_PARAMETER, 'amber/AMBER_rotamers.dat') )                            #
try:                                                                                                                                #
    os.remove(TRAJECTORY)                                                                                                           #
except:                                                                                                                             #
    pass                                                                                                                            #
#-----------------------------------------------------------------------------------------------------------------------------------#
system.energy(log = True)

'''
pdbs          = [
                 os.path.join(PEPDICE_EXAMPLES ,'data/alpha/1GAB/1gab_0001.pdb'),
                 os.path.join(PEPDICE_EXAMPLES ,'data/alpha/1GAB/1gab_0002.pdb'),
                 os.path.join(PEPDICE_EXAMPLES ,'data/alpha/1GAB/1gab_0003.pdb'),
                 os.path.join(PEPDICE_EXAMPLES ,'data/alpha/1GAB/1gab_0004.pdb'),
                 os.path.join(PEPDICE_EXAMPLES ,'data/alpha/1GAB/1gab_0005.pdb'),
                 os.path.join(PEPDICE_EXAMPLES ,'data/alpha/1GAB/1gab_0006.pdb'),
                 os.path.join(PEPDICE_EXAMPLES ,'data/alpha/1GAB/1gab_0007.pdb'),
                 os.path.join(PEPDICE_EXAMPLES ,'data/alpha/1GAB/1gab_0008.pdb'),
                 os.path.join(PEPDICE_EXAMPLES ,'data/alpha/1GAB/1gab_0009.pdb'),
                 os.path.join(PEPDICE_EXAMPLES ,'data/alpha/1GAB/1gab_0010.pdb'),
                 os.path.join(PEPDICE_EXAMPLES ,'data/alpha/1GAB/1gab_0011.pdb'),
                 os.path.join(PEPDICE_EXAMPLES ,'data/alpha/1GAB/1gab_0012.pdb'),
                 os.path.join(PEPDICE_EXAMPLES ,'data/alpha/1GAB/1gab_0013.pdb'),
                 os.path.join(PEPDICE_EXAMPLES ,'data/alpha/1GAB/1gab_0014.pdb'),
                 os.path.join(PEPDICE_EXAMPLES ,'data/alpha/1GAB/1gab_0015.pdb'),
                 os.path.join(PEPDICE_EXAMPLES ,'data/alpha/1GAB/1gab_0016.pdb'),
                 os.path.join(PEPDICE_EXAMPLES ,'data/alpha/1GAB/1gab_0017.pdb'),
                 os.path.join(PEPDICE_EXAMPLES ,'data/alpha/1GAB/1gab_0018.pdb'),
                 os.path.join(PEPDICE_EXAMPLES ,'data/alpha/1GAB/1gab_0019.pdb'),
                 os.path.join(PEPDICE_EXAMPLES ,'data/alpha/1GAB/1gab_0020.pdb')]
'''



pdbs = [
        os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/native.pdb'),   
        os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy0_1.pdb'),   
        os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy1_1.pdb'),   
        os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy1_7.pdb'),   
        os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy1_8.pdb'),   
        os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy1_3.pdb'),   
        os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy1_4.pdb'),   
        os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy1_15.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy1_19.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy1_9.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy1_10.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy1_11.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy1_12.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy1_13.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy1_2.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy1_5.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy1_6.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy1_18.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy1_14.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy1_16.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy1_17.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy1_41.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy1_42.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy1_57.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy1_59.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy1_21.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy1_23.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy1_24.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy1_35.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy1_22.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy1_27.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy1_31.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy1_33.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy1_34.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy1_36.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy1_37.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy1_39.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy1_40.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy1_43.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy1_48.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy1_49.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy1_53.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy1_62.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy1_64.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy1_65.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy2_12.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy2_237.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy2_240.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy2_241.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy2_246.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy2_253.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy2_254.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy2_232.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy2_259.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy2_316.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy2_303.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy2_300.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy3_1.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy3_129.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy3_131.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy3_158.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy3_162.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy3_164.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy3_166.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy3_140.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy3_144.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy3_155.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy3_157.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy3_176.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy3_179.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy3_183.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy3_186.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy3_187.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy3_193.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy3_238.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy3_235.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy3_240.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy3_242.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy3_276.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy3_272.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy3_275.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy3_282.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy3_311.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy3_310.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy4_5.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy4_108.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy4_130.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy4_136.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy4_158.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy4_176.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy4_197.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy4_200.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy4_206.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy4_228.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy4_243.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy4_229.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy4_235.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy4_236.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy4_252.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy4_226.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy4_234.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy4_238.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy5_19.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy5_50.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy5_107.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy5_119.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy5_144.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy6_10.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy6_276.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy6_297.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy6_298.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy7_2.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy7_59.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy7_72.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy7_134.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy7_124.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy7_126.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy7_121.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy7_135.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy7_173.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy7_175.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy8_1.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy9_40.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy9_170.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy9_163.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy9_173.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy9_179.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy9_200.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy9_211.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy9_209.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy9_214.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy9_207.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy9_234.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy9_241.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy9_243.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy9_250.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy9_269.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy9_280.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy9_281.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy9_273.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy9_277.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy9_270.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy9_271.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy9_275.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy9_296.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy9_303.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy10_4.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy10_86.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy10_182.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy10_222.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy10_272.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy10_279.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy11_4.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy12_5.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy12_6.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy12_8.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy12_35.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy12_118.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy13_1.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy14_14.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy15_1.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy16_1.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy17_18.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy18_9.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy19_5.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy19_172.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy20_1.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy21_1.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy21_147.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy22_18.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy23_17.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy24_22.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy25_3.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy25_68.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy25_162.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy26_6.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy26_53.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy26_55.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy26_64.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy27_9.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy27_10.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy27_19.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy27_15.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy27_4.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy27_8.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy27_22.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy27_24.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy27_44.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy27_68.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy27_78.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy27_80.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy28_6.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy28_8.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy28_11.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy28_56.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy28_64.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy29_6.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy29_89.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy29_94.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy29_104.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy29_80.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy29_81.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy29_88.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy29_95.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy29_96.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy30_1.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy31_17.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy32_16.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy32_98.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy33_5.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy33_140.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy33_145.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy33_151.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy33_143.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy33_154.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy33_157.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy34_1.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy34_37.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy34_56.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy34_98.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy34_107.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy34_113.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy34_88.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy34_145.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy34_205.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy34_193.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy34_199.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy34_223.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy34_226.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy34_194.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy34_195.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy34_202.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy34_209.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy34_212.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy34_238.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy35_11.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy36_1.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy37_22.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy38_1.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy39_3.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy40_5.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy41_5.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy41_41.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy41_30.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy41_56.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy41_76.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy41_82.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy41_134.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy42_8.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy42_50.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy43_2.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy44_3.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy44_88.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy45_1.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy46_2.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy47_6.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy47_91.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy47_95.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy48_29.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy49_7.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy49_205.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy49_220.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy49_222.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy49_235.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy50_3.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy51_26.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy52_2.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy52_68.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy53_37.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy53_159.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy53_164.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy53_166.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy54_3.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy55_22.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy56_25.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy57_29.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy58_6.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy59_15.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy59_79.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy59_62.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy60_47.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy61_7.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy62_1.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy62_106.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy63_3.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy64_2.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy64_144.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy65_25.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy66_2.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy66_85.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy66_117.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy66_121.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy66_125.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy66_127.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy66_137.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy66_171.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy66_193.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy66_201.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy66_202.pdb'),   
        #os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy67_14.pdb'),   
        os.path.join(PEPDICE_EXAMPLES, 'data/alpha/1GAB/1GAB/decoy68_18.pdb'),]


system  = build_fragment_library_from_pdbs (
                                            molecule             = system ,
                                            frag_size            = 52     ,
                                            number_of_fragments  = 1      ,
                                            pdblist              = pdbs   ,
                                            )

import pickle

pickle.dump(system.fragments, open( "1gab_fragments.p", "wb" ) )


