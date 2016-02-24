#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  MDLigRefiner_tleap.py
#  
#  Copyright 2015 farminf <farminf@farminf-3>
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



import sys
import os
import pprint as pprint
import subprocess

'''
class ANTECHAMBER:
    """ Class doc """

    def __init__ (self                , 
                 filein        = None ,
                 charge        = 0    , 
                 multiplicity  = 1    , 
                 charge_method = 'bcc', 
                 ):
                     
        self.filein        = filein.strip()
                           
        file_in2           = filein.split("/")
        file_in2           = file_in2[-1]
        file_in2           = file_in2.split(".")  
                           
        self.file_type     = file_in2[-1]
        self.file_name     = file_in2[-2]
        self.charge        = str(charge      )
        self.multiplicity  = str(multiplicity)
        self.charge_method = charge_method

        self.mol2file      = None

    def ANTECHAMBER_run_antechamber(self):

        AMBERTOOLS_outputs = "AMBERTOOLS_outputs"

        if not os.path.exists ( AMBERTOOLS_outputs ): os.mkdir ( AMBERTOOLS_outputs )
        try:
            os.rename(AMBERTOOLS_outputs+"/"+new_ligand_name+".mol2",AMBERTOOLS_outputs+"/"+new_ligand_name+".mol2.backup")
        except:
            a = None

        print AMBERTOOLS_outputs
        print self.file_name 
        print self.file_type 
        print self.charge_method
        print self.filein
        print "Charge model: ",self.charge_method
        print "Charge:       ",self.charge
        print "Multiplicity: ",self.multiplicity

         
        command  =  'antechamber '
        command +=  " -i  " + self.filein    + " "
        command +=  " -fi " + self.file_type + " "

        command +=  " -o  " + AMBERTOOLS_outputs + "/" + self.file_name + ".mol2 "+" -fo mol2 "
        command +=  " -rn " + self.file_name     + " "
        command +=  " -c  " + self.charge_method + " "
        command +=  " -nc " + str(self.charge) + " -m " + str(self.multiplicity)
        print command
        os.system (command)

        #os.system("antechamber -i "+AMBERTOOLS_outputs+"/"+new_ligand_name+".pdb -fi pdb -o "+AMBERTOOLS_outputs+"/"+new_ligand_name+".mol2 "+" -fo mol2  -c "+ charge_method + " -nc "+str(charge)+" -m "+ str(multiplicity))

        print "done"							
                
        try:
            os.rename("ANTECHAMBER_AC.AC0", AMBERTOOLS_outputs+ "/ANTECHAMBER_AC.AC0")
        except:
            a = None
        try:
            os.rename("ANTECHAMBER_AM1BCC.AC", AMBERTOOLS_outputs+ "/ANTECHAMBER_AM1BCC.AC")
        except:
            a = None
        try:
            os.rename("ANTECHAMBER_AC.AC", AMBERTOOLS_outputs+ "/ANTECHAMBER_AC.AC")
        except:
            a = None
        try:
            os.rename("ANTECHAMBER_AM1BCC_PRE.AC", AMBERTOOLS_outputs+ "/ANTECHAMBER_AM1BCC_PRE.AC")
        except:
            a = None
        try:
            os.rename("ANTECHAMBER_BOND_TYPE.AC", AMBERTOOLS_outputs+ "/ANTECHAMBER_BOND_TYPE.AC")
        except:
            a = None
        try:
            os.rename("ANTECHAMBER_BOND_TYPE.AC0", AMBERTOOLS_outputs+ "/ANTECHAMBER_BOND_TYPE.AC0")
        except:
            a = None
        try:
            os.rename("ATOMTYPE.INF", AMBERTOOLS_outputs+ "/ATOMTYPE.INF")
        except:
            a = None
        try:
            os.rename("sqm.in", AMBERTOOLS_outputs+ "/sqm.in")
        except:
            a = None
        try:
            os.rename("sqm.out", AMBERTOOLS_outputs+ "/sqm.out")
        except:
            a = None

         
        print "now runing parmchk\n"
        os.system("parmchk -i " + AMBERTOOLS_outputs + "/" + self.file_name + ".mol2 " + " -f mol2 -o "+AMBERTOOLS_outputs+"/"+ self.file_name +".frcmod")
        file_in = open(AMBERTOOLS_outputs+"/"+self.file_name+".frcmod", 'r')

        for line in file_in:
            print line
        print "\n done \n"
'''

class TLEAP:
    """ Class doc """
    def __init__ (self, WORKSPACE = None):
        """ Class initialiser """

        #    - WORKSPACE - 
        if WORKSPACE == None:
            self.WORKSPACE  =  os.path.join(os.getcwd(), 'workspace')
        else:
            self.WORKSPACE  =  WORKSPACE

        try:
            if not os.path.isdir(self.WORKSPACE):
                os.mkdir(self.WORKSPACE)
        except:
            print WORKSPACE, ' already exist'
            
            
        self.AMBERHOME = os.environ.get('AMBERHOME')

        print 'AMBERHOME found: ', self.AMBERHOME


        oldffdir = os.path.join(self.AMBERHOME, 'dat/leap/cmd/oldff')
        self.old_ff_list = os.listdir(oldffdir)

        newffdir  = os.path.join(self.AMBERHOME, 'dat/leap/cmd')
        self.new_ff_list  = os.listdir(newffdir)

        self.force_fields = {
                       #
                       } 

        for ff in self.new_ff_list:
            self.force_fields[ff] = os.path.join(self.AMBERHOME, 'dat/leap/cmd', ff)
        for ff in self.old_ff_list:
            self.force_fields[ff] = os.path.join(self.AMBERHOME, 'dat/leap/cmd/oldff', ff)

        self.lib_list    = []
        self.prep_list   = []
        self.mol2_list   = []
        self.frcmod_list = []

        for i in self.force_fields:
            print i, self.force_fields[i]

    def TLEAP_script_maker (self                      , 
                            pdbin         = None      , # pdb fullpath (including modifications) str	
                            ligand_list   = []        , # mol2, prep, lib and frcmod files
                            leaprc        = 'leaprc'  , 
                            system_name   = 'mysystem', # output files           
                            ff_model      = None      , # ff model - ff99SB                      str
                            ff_model_gly  = None      , # glycan type  ex glycam 06              str
                            ff_model_gaff = None      , # include gaff?                       True/False
                            water_model   = None      , # water - TIP3P                          str
                            water_box     = None      , # include waterbox ?                  True/False
                            add_ions      = None      , # include ions?                       True/False
                            neutralize    = None      , # neutralize total charge ?           True/False  
                            positive_type = None      , # add ions Ex Na+                        str
                            negative_type = None      , # add ions Ex Cl-                        str
                            positive_num  = None      , # number of positive ions                int
                            negative_num  = None      ,
                            link_atoms    = []        ):                                     
        """ Function doc """

        #---------------------------------#
        #           force_field           #
        #---------------------------------#
        force_field = self.force_fields[ff_model]


        text=[]
        # ff_model

        text.append("source "+force_field+"\n")
        print "source       "+force_field+"\n"		



        #---------------------------------#
        #       glycan_force_field        #
        #---------------------------------#		
        if ff_model_gly != None:
            glycan_force_field = self.force_fields[ff_model_gly]
            #print "source  "+glycan_force_field+"\n"
            text.append("\nsource  "+glycan_force_field+"\n")



        #---------------------------------#
        #       glycan_force_field        #
        #---------------------------------#		
        if ff_model_gaff == True:
            text.append("source leaprc.gaff\n")
            #print "source leaprc.gaff\n"
            

        #------------------------------------------------------------------#
        #            ligand_list mol2, lib, prep and frcmod files          #
        #------------------------------------------------------------------#	
        for lig in ligand_list:
            print lig
            lig      = lig.strip()
            file_in  = lig.strip()
            file_in2 = lig.split("/")
            file_in2 = file_in2[-1]
            file_in2 = file_in2.split(".")
            #print "residue name: ", file_in2[-2]
            #print "file type: "   , file_in2[-1]
            file_type = file_in2[-1]
            print file_type
            if file_type == "mol2":
                text.append(file_in2[-2]+" = loadmol2 "+file_in+"\n")
                print "mol2 aqui"
            
            if file_type  == "frcmod":
                text.append("loadamberparams "+file_in+"\n")
                print "mol2 aqui"
            
            if file_type  == "prep":
                text.append("loadamberprep "+file_in+"\n")
                print "mol2 aqui"
            
            if file_type  == "lib":
                text.append("loadoff "+file_in+"\n")
                print "mol2 aqui"
            
            #print 'fim'
        #----------------------#
        #        Ions FF       #
        #----------------------#                  
                    
        if ff_model in self.new_ff_list:
            #text.append("\nlist\n")
            text.append("\nloadAmberParams frcmod.ionsjc_tip3p\n")


        print 'loadpdb'
        #----------------------#
        #        loadpdb       #
        #----------------------#                  
        text.append("system = loadpdb " + pdbin +"\n")
        #text.append("system = loadpdb "+dataPath+"/file_teste.pdb\n")

        for link in link_atoms:
            text.append(link + "\n")


        #water_model = App.get_object('cbox_tleap_solvate_water_model').get_active_text()
        if water_model != None:
            text.append("solvatebox system " + water_model +"BOX " + water_box + "\n")

        if add_ions == True:	
            if neutralize == True:
                text.append("addions system Cl- 0\n")
                text.append("addions system Na+ 0\n")
            else:
                #print "use 'just neutralize charges' "
                text.append("addions system " + negative_type + " " + negative_num + " \n")
                text.append("addions system " + positive_type + " " + positive_num + " \n")


                #positive_type = add ions Ex Na+                str
                #negative_type = add ions Ex Cl-                str
                #positive_num  = number of positive ions        int
                #negative_num): 

        #text.append("saveamberparm system " +data_path+"/"+system_name+".top " +data_path+"/"+system_name+".crd \n")
        text.append("saveamberparm system " +system_name+".top " +system_name+".crd \n")

        text.append("quit \n")
        arq = open(leaprc, 'w')
        arq.writelines(text)
        arq.close()    

    def TLEAP_run (self                     ,
                   inputfile  = 'leaprc'    ,
                   outputfile = 'tleap.log'
                   ):
        """ Function doc """
        
        null_file = open(os.devnull, "w")
        commandTleap = 'tleap -f '+inputfile+' > '+outputfile
        subprocess.call(commandTleap.split(), stdout=null_file, stderr=null_file)
        
        #os.system('tleap -f '+inputfile+' > '+outputfile)

        def TLEAP_logfile_reader (self):
            """ Function doc """
            pass

    def TLEAP_convert_topology_from_amber12_to_amber11 (self, filein = None, fileout = None):
        filein = open(filein, 'r')
        text   = []
        print_line = True

        for line in filein:
            line2 = line.split()
            try:
                if line2[0] == '%FLAG':
                    if   line2[1] == 'ATOMIC_NUMBER':
                        print 'excluding flag:', line/home/farminf/pDynamoWorkSpace/Ramon_Nov_23_2015/MOPAC_files/tmp/system0.arc
                        print_line = False

                    elif   line2[1] == 'SCEE_SCALE_FACTOR':
                        print 'excluding flag:', line
                        print_line = False

                    elif   line2[1] == "SCNB_SCALE_FACTOR":
                        print 'excluding flag:', line
                        print_line = False			

                    elif   line2[1] == 'IPOL':
                        print 'excluding flag:', line
                        print_line = False
            
                    else:
                        print_line = True	
                        #print print_line
            except:
                a= None
            if print_line == True:
                text.append(line)

        fileout = open(fileout, 'w')
        fileout.writelines(text)
        fileout.close()

        ligand_list = [ "/home/fernando/Dropbox/Python/source_files/2AX.frcmod   ",
                "/home/fernando/Dropbox/Python/source_files/2AX.mol2     ",
                "/home/fernando/Dropbox/Python/source_files/2AX.pdbst    ",
                "/home/fernando/Dropbox/Python/source_files/nad.frcmod   ",
                "/home/fernando/Dropbox/Python/source_files/nadh.prep    ",
                "/home/fernando/Dropbox/Python/source_files/nadp+.prep   ",
                "/home/fernando/Dropbox/Python/source_files/nadph.frcmod ",
                "/home/fernando/Dropbox/Python/source_files/nadph.prep   "]

        pdbin = '/home/fernando/Dropbox/Python/source_files/receptor_CW_2AX.pdb'




class ANTECHAMBER:
    """ Class doc """

    def __init__ (self               ,
                  i  = None          ,  #   input file name
                  fi = None          ,  #   input file format
                  o  = None          ,  #   output file name
                  fo = None          ,  #   output file format
                  c  = None          ,  #   charge method
                  cf = None          ,  #   charge file name
                  nc = None          ,  #   net molecular charge (int)
                  a  = None          ,  #   additional file name
                  fa = None          ,  #   additional file format
                  ao = None          ,  #   additional file operation
                                        #       crd   : only read in coordinate
                                        #       crg   : only read in charge
                                        #       radius: only read in radius
                                        #       name  : only read in atom name
                                        #       type  : only read in atom type
                                        #       bond  : only read in bond type 
                                    
                  m  = None          ,  #   multiplicity (2S+1), default is 1
                  rn = 'MOL'          ,  #       residue name, overrides input file, default is MOL
                                    
                  rf = None          ,  #   residue toplogy file name in prep input file,
                                        #       default is molecule.res
                                    
                  ch = None          ,  #   check file name for gaussian, default is 'molecule'
                  ek = None          ,  #   mopac or sqm keyword, inside a pair of quotes
                  gk = None          ,  #   gaussian job keyword, inside a pair of quotes
                  gm = None          ,  #   gaussian memory keyword, inside a pair of quotes, such as"%mem=1000MB"
                  gn = None          ,  #   gaussian number of processors keyword, inside a pair of quotes, such as "%nproc=8"
                  gv = None          ,  #   add keyword to generate gesp file (for Gaussian 09 only)
                                        #       1    : yes
                                        #       0    : no, the default
                                    
                  ge = None          ,  #   gaussian esp file generated by iop(6/50=1), default is g09.gesp
                  df = None          ,  #   am1-bcc precharge flag, 2 - use sqm(default); 0 - use mopac
                  at = None          ,  #   atom type, can be gaff (the default), amber (for PARM94/99/99SB), bcc and sybyl
                  du = None          ,  #   fix duplicate atom names: yes(y)[default] or no(n)
                  bk = None          ,  #   4-character component Id, for ccif
                  an = None          ,  #   adjust atom names: yes(y) or no(n)
                                        #   the default is 'y' for 'mol2' and 'ac' and 'n' for the other formats 
                                    
                  j  = None          ,  #   atom type and bond type prediction index, default is 4 
                                        #       0    : no assignment
                                        #       1    : atom type 
                                        #       2    : full  bond types 
                                        #       3    : part  bond types 
                                        #       4    : atom and full bond type 
                                        #       5    : atom and part bond type 
                                    
                  s  = None          ,  #   status information: 0(brief), 1(default) or 2(verbose)
                  eq = None          ,  #   equalizing atomic charge, default is 1 for '-c resp' and '-c bcc' and 0 for the other charge methods 
                                        #       0    : no use
                                        #       1    : by atomic paths 
                                        #       2    : by atomic paths and structural information, i.e. E/Z configurations 
                                    
                  pf = None          ,  #   remove intermediate files: yes(y) or no(n)[default]
                  pl = None          ,  #   maximum path length to determin equivalence of atomic charges for resp and bcc,
                                        #       the smaller the value, the faster the algorithm, default is -1 (use full length),
                                        #       set this parameter to 10 to 30 if your molecule is big (# atoms >= 100)
                                    
                  directory = None   ): #   output directory 

                  #Use 'antechamber -L' to list the supported file formats and charge methods

        self.i  = i  
        self.fi = fi 
        self.o  = o  
        self.fo = fo 
        self.c  = c  
        self.cf = cf 
        self.nc = nc 
        self.a  = a  
        self.fa = fa 
        self.ao = ao 
        self.m  = m  
        self.rn = rn 
        self.rf = rf 
        self.ch = ch 
        self.ek = ek 
        self.gk = gk 
        self.gm = gm 
        self.gn = gn 
        self.gv = gv 
        self.ge = ge 
        self.df = df 
        self.at = at 
        self.du = du 
        self.bk = bk 
        self.an = an 
        self.j  = j  
        self.s  = s  
        self.eq = eq 
        self.pf = pf 
        self.pl = pl 
        self.i  = i  
        
        if directory == None:
            self.AMBERTOOLS_outputs = os.getcwd()
        else:
            self.AMBERTOOLS_outputs = directory
            
            
    def build_AMBERTOOLS_outputs_diretory(self, AMBERTOOLS_outputs = "AMBERTOOLS_outputs"):
        """ Function doc """
        #AMBERTOOLS_outputs = "AMBERTOOLS_outputs"
        if not os.path.exists ( AMBERTOOLS_outputs ): os.mkdir ( AMBERTOOLS_outputs )
        self.AMBERTOOLS_outputs = AMBERTOOLS_outputs
        
    def change_AMBERTOOLS_outputs_directory (self, new_directory = None):
        """ Function doc """
        if new_directory:
            self.AMBERTOOLS_outputs = new_directory
    
    def backup_AMBERTOOLS_outputs (self):
        """ Function doc """
        try:
            os.rename("ANTECHAMBER_AC.AC0", self.AMBERTOOLS_outputs+ "/ANTECHAMBER_AC.AC0")
        except:
            a = None
        try:
            os.rename("ANTECHAMBER_AM1BCC.AC", self.AMBERTOOLS_outputs+ "/ANTECHAMBER_AM1BCC.AC")
        except:
            a = None
        try:
            os.rename("ANTECHAMBER_AC.AC", self.AMBERTOOLS_outputs+ "/ANTECHAMBER_AC.AC")
        except:
            a = None
        try:
            os.rename("ANTECHAMBER_AM1BCC_PRE.AC", self.AMBERTOOLS_outputs+ "/ANTECHAMBER_AM1BCC_PRE.AC")
        except:
            a = None
        try:
            os.rename("ANTECHAMBER_BOND_TYPE.AC", self.AMBERTOOLS_outputs+ "/ANTECHAMBER_BOND_TYPE.AC")
        except:
            a = None
        try:
            os.rename("ANTECHAMBER_BOND_TYPE.AC0", self.AMBERTOOLS_outputs+ "/ANTECHAMBER_BOND_TYPE.AC0")
        except:
            a = None
        try:
            os.rename("ATOMTYPE.INF", self.AMBERTOOLS_outputs+ "/ATOMTYPE.INF")
        except:
            a = None
        try:
            os.rename("sqm.in", self.AMBERTOOLS_outputs+ "/sqm.in")
        except:
            a = None
        try:
            os.rename("sqm.out", self.AMBERTOOLS_outputs+ "/sqm.out")
        except:
            a = None

    def run_antechamber(self):
        print self.AMBERTOOLS_outputs
        print self.i
        print self.fi
        print "Charge model: ",self.c
        print "Charge:       ",self.nc
        print "Multiplicity: ",self.m

        command  =  'antechamber '
        command +=  " -i  " + self.i  + " "
        command +=  " -fi " + self.fi + " "
        command +=  " -o  " + self.o  + " "
        command +=  " -fo " + self.fo + " "
        command +=  " -rn " + self.rn + " "
        command +=  " -c  " + self.c  + " "
        command +=  " -nc " + str(self.nc) + " " + " -m " + str(self.m)

        print command
        os.system (command)

        #os.system("antechamber -i "+AMBERTOOLS_outputs+"/"+new_ligand_name+".pdb -fi pdb -o "+AMBERTOOLS_outputs+"/"+new_ligand_name+".mol2 "+" -fo mol2  -c "+ charge_method + " -nc "+str(charge)+" -m "+ str(multiplicity))

        print "done"							
                
    def run_parmchk (self, 
                     i = None, 
                     f = 'mol2', 
                     o = 'output.frcmod'):
        """ Function doc """
        print "now runing parmchk\n"
        
        try:
            if i:
                command  =  "parmchk " 
                command +=  " -i " + i
                command +=  " -f " + f
                command +=  " -o " + o

                os.system (command)
                file_in = open(o, 'r')
                for line in file_in:
                    print line
                print "\n done \n"
        except:
            print 'fail: ', command
            


def main():
    
    source = '/home/farminf/Dropbox/Python/source_files/'
    lista = os.listdir(source)
    
    
    #for File in lista:
    #    
    #    Type =  File.split('.')
    #    name =  Type[0]
    #    if Type[-1] == 'pdb':
    #        
    #
    #        a = ANTECHAMBER(
    #                        i  = os.path.join(source,File)          ,
    #                        fi = Type[-1]                           ,
    #                        
    #                        o  = os.path.join(source,name + '.mol2'),
    #                        fo = 'mol2'                             ,
    #                        rn =  name                              ,
    #                        c  = 'bcc'                              ,
    #                        nc = '0'                                ,
    #                        m  = '1'                                ,
    #                        )
    #        
    #        a.run_antechamber()
    #        a.run_parmchk  (i = os.path.join(source,name + '.mol2'), 
    #                        f = 'mol2'               , 
    #                        o = os.path.join(source,name + '.frcmod'))
    #        
            
    pdbin         =  '/home/farminf/Dropbox/Python/source_files/receptor_CW_2AX.pdb'
    ligand_list   = ['/home/farminf/Dropbox/Python/source_files/nad.frcmod'  ,
                     '/home/farminf/Dropbox/Python/source_files/nadh.prep'   ,
                     '/home/farminf/Dropbox/Python/source_files/2AX.frcmod'  ,
                     '/home/farminf/Dropbox/Python/source_files/2AX.mol2'    ,
                     ]
    
    
    system_name = 'mysystem'
    a = TLEAP()
    a.TLEAP_script_maker (
                          pdbin         = pdbin                ,
                          ligand_list   = ligand_list          ,
                          system_name   = system_name          ,
                          ff_model      = 'leaprc.ff99SB'      ,
                          ff_model_gly  = 'leaprc.GLYCAM_06EP' ,
                          ff_model_gaff = True      ,
                          water_model   = 'TIP3P'   ,
                          water_box     = '10'      ,
                          add_ions      = True      ,
                          neutralize    = True      ,
                          positive_type = 'Na+'     ,
                          negative_type = 'Cl-'     ,
                          positive_num  = None      ,
                          negative_num  = None      ,
                          link_atoms    = []        )
    a.TLEAP_run()
    a.TLEAP_convert_topology_from_amber12_to_amber11 ( filein = system_name+'.top' , fileout = system_name+'11.top')
    return 0

if __name__ == '__main__':
	main()

