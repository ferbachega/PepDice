from ParseNaMDLog import ParseNaMDLog
from ParseGMXLog  import ParseGMXLog
from ParseAMBERLog import ParseAMBERLog
from CRDFiles      import save_CRD_to_file
from Geometry    import  distance_ab

import os
import sys



def compute_vdw_ij (atom_i, atom_j):
    """ Function doc """ #-23085572255.9
    A    = 4
    B    = 1
    R_ab = distance_ab (atom_i, atom_j)
    E_ab = A*((R_ab**-12) - B*(atom_i.epsilon*atom_j.epsilon)*(R_ab**-6))
    return E_ab



def compute_AB_vdw_ij (atom_i, atom_j):
    """ Function doc """ #-23085572255.9
    
    A    = 1
    B    = 1
    
    if atom_i.AB + atom_j.AB == 'AA':   # apolares
        C = 1
    elif atom_i.AB + atom_j.AB == 'BB': # polares
        C = 0.5
    else:
        C = -0.5
    
    R_ab = distance_ab (atom_i, atom_j)
    E_ab = A*((R_ab**-12) - C*(R_ab**-6))
    return E_ab

def compute_AB_energy (molecule = None):
    """ Function doc """
    total_E = 0
    #----------------------- atom i -----------------------------#
    for residue_i in molecule.residues:
        for atom_i in residue_i.atoms:
            
            if atom_i.name == 'CA':
                if residue_i.name in ['ALA','ILE','LEU','MET','VAL','PRO','GLY','CYS']:
                    atom_i.AB = 'A'
                else:
                    atom_i.AB = 'B'
                #------------------- atom j ---------------------#
                for residue_j in molecule.residues:
                    if residue_j.id == residue_i:
                        pass
                    else:
                        
                        for atom_j in  residue_j.atoms:
                            
                            if atom_i.id == atom_j.id:
                                pass
                            
                            else:
                                if atom_j.name == 'CA':
                                    if residue_j.name in ['ALA','ILE','LEU','MET','VAL','PRO','GLY','CYS']:
                                        atom_j.AB = 'A'
                                    else:
                                        atom_j.AB = 'B'
                                    vdw =  compute_AB_vdw_ij (atom_i, atom_j)
                                    total_E += vdw
    return total_E





class Energy:
    """ Class doc """
    
    def __init__ (self):
        """ Class initialiser """
        pass

    def energy(self, 
               log                       = False, 
               pn                        = 1    ,  #process number # used in multiprocess 
               external_coordinates      = False, 
               external_coordinates_type = 'pdb',
               external_coordinates_file = None , 
               #energy_AB                 = True ,
               ):
                    
        """ Function doc """
        bond     = self.bond    
        angle    = self.angle   
        dihed    = self.dihed   
        imprp    = self.imprp   
        elect    = self.elect   
        vdw      = self.vdw     
        boundary = self.boundary
        esurf    = self.esurf   
        egb      = self.egb     
        

        if  self.ff_type == 'charmm':
            write_NaMD_input_file (molecule = self, Type='energy', pn = pn)
            save_PDB_to_file      (molecule = self, filename='SinglePoint'+str(pn)+'.pdb')
       
            os.system('namd2 SinglePoint.namd > SinglePoint'+str(pn)+'.log')
            
            BOND , ANGLE , DIHED , IMPRP , ELECT , VDW , BOUNDARY = ParseNaMDLog('SinglePoint'+str(pn)+'.log', log=log)
            energy = (BOND*bond + ANGLE*angle + DIHED*dihed + 
                      IMPRP*imprp + ELECT*elect + VDW*vdw + BOUNDARY*boundary)
            
            return energy


        if self.ff_type == 'gmx':
            import subprocess
            
            pn  = pn
            if external_coordinates:
                command1   = 'mdrun -s '+ self.tpr+' -rerun '+ external_coordinates_file +' -e ener_'+str(pn)
                
                #os.system('mdrun -s '+ self.tpr+' -rerun '+ external_coordinates_file +' -e ener_'+str(pn))
                #mdrun -s test.tpr -rerun complex_2.pdb
                command2   = 'echo 1 2 3 4 5 8 33 | g_energy -f ener_'+str(pn)+'.edr -s '+ self.tpr + ' -o SP_'+str(pn)
                #print command2
                
                null_file = open(os.devnull, 'w')
                subprocess.call(command1.split(), stdout = null_file, stderr = null_file)
                os.system('echo 1 2 3 4 5 6 10 | g_energy -f ener_'+str(pn)+'.edr -s '+ self.tpr + ' -o SP_'+str(pn))
                #subprocess.call(command2.split(), stdout = null_file, stderr = null_file)
                energy_list = ParseGMXLog('SP_'+str(pn)+'.xvg', log=log)
                
                os.system('rm *#')
                return energy_list['Potential']

        
        if  self.ff_type == 'amber':
            # transformar numa funcao
            pn  = pn
            
            write_AMBER_input_file(molecule = self, Type='energy', pn= pn)
            save_CRD_to_file      (molecule = self, filename='SinglePoint'+str(pn)+'.crd')
            
            os.system('sander -O -i SinglePoint'+str(pn)+'.in -c SinglePoint'+str(pn)+'.crd -o SinglePoint'+str(pn)+'.log -p ' + self.top)
            
            #sander -O -i energy.in  -o   energy.log -p 7tim.top -c 7tim.crd
            
            energy_list = ParseAMBERLog('SinglePoint'+str(pn)+'.log', log=log)
            #print BOND , ANGLE , DIHED , IMPRP , ELECT , VDW , BOUNDARY
            
            energy = 0 
            energy += energy_list["ESURF"]       * esurf
            #energy += energy_list["RESTRAINT"]
            energy += energy_list["EGB"]         * egb
            energy += energy_list["EELEC"]       * elect

            if energy_list["VDWAALS"] == None:
                #energy_list["VDWAALS"] = 99999999999999999999999
                #energy += energy_list["VDWAALS"] * vdw
                return None
            else:
                energy += energy_list["VDWAALS"] * vdw
            
            #energy += energy_list["EEL"]
            #energy += energy_list["NB"]
            energy += energy_list["DIHED"] * dihed
            energy += energy_list["ANGLE"] * angle
            energy += energy_list["BOND"]  * bond
            
            if self.AB != 0.0:
                energy +=   compute_AB_energy(molecule = self)*self.AB
            
            #print 'energy_AB:', energy_AB 
            
            if log:
                from pprint import pprint
                pprint(energy_list)
            return energy


        if  self.ff_type == 'Calpha_model':
            total_E = 0
            #----------------------- atom i -----------------------------#
            for residue_i in self.residues:
                for atom_i in residue_i.atoms:
                    if atom_i.name == 'CA':
                        #------------------- atom j ---------------------#
                        for residue_j in self.residues:
                            if residue_j.id == residue_i:
                                pass
                            else:
                                
                                for atom_j in  residue_j.atoms:
                                    
                                    if atom_i.id == atom_j.id:
                                        pass
                                    
                                    else:
                                        if atom_j.name == 'CA':
                                            vdw =  compute_vdw_ij (atom_i, atom_j)
                                            #print atom_i.id, residue_i.name, atom_i.name,atom_i.epsilon, atom_j.id, residue_j.name, atom_j.name, atom_j.epsilon, distance_ab (atom_i, atom_j), vdw
                                            total_E += vdw

            if log == True:
                print 'VdW energy: ',total_E
            return total_E




def write_AMBER_input_file (molecule=None, Type='energy', pn = 1):
    text = """  compute single-point energy 
 &cntrl
   cut=999.0, igb=1, saltcon=0.2, gbsa=1, rgbmax = 999.00000, surften = 0.010,
   ntpr=1,
   nstlim = 0, dt=0.002,
   ntt=1, tempi=300.0, temp0=300.0, tautp=2.0,
   ntx=1, irest=0, ntb=0,
 &end
eof"""
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
            
        #if _type == 'dynamics':
        #text += 'minimize 100                                                             \n'    
       #text += 'minimize            1000                                                 \n'

    #output_file = open('SinglePoint.namd', "w")
    #output_file.write(text)
    #output_file.close()
