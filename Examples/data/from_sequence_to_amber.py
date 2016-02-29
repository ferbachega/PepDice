#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  from_seq2struc.py
#  
#  Copyright 2014 Labio <labio@labio-XPS-8300>
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
import random
import os

one_letter   ={'VAL':'V', 'ILE':'I', 'LEU':'L', 'GLU':'E', 'GLN':'Q', \
               'ASP':'D', 'ASN':'N', 'HIS':'H', 'TRP':'W', 'PHE':'F', 'TYR':'Y',    \
               'ARG':'R', 'LYS':'K', 'SER':'S', 'THR':'T', 'MET':'M', 'ALA':'A',    \
               'GLY':'G', 'PRO':'P', 'CYS':'C'}

three_letter ={'V':'VAL', 'I':'ILE', 'L':'LEU', 'E':'GLU', 'Q':'GLN',               \
			   'D':'ASP', 'N':'ASN', 'H':'HIS', 'W':'TRP', 'F':'PHE', 'Y':'TYR',    \
			   'R':'ARG', 'K':'LYS', 'S':'SER', 'T':'THR', 'M':'MET', 'A':'ALA',    \
			   'G':'GLY', 'P':'PRO', 'C':'CYS'}

#sequence              = 'MHVTQSSSAITPGQTAELYPGDIKSVLLTAEQIQARIAELGEQIGNDYRELSATTGQ'
#ss_sequence           = 'CCEEECCCCCCCCCCCCCCHHHHHHEECCHHHHHHHHHHHHHHHHHHHHHCCCCCCC'


#-------------------------------------------------------------------------------#
#                     pdb 1ZDD      - double helix hairpin                      #
#-------------------------------------------------------------------------------#
sequence              = 'FNMQCQRRFYEALHDPNLNEEQRNAKIKSIRDDCG'
ss_sequence           = 'CCHHHHHHHHHHHCCCCCCHHHHHHHHHHHHCCCC'

sequence              = 'EQYTAKYKGRTFRNEKELRDFIEKFKGR'
ss_sequence           = 'CCCCCCCCCCCCCCHHHHHHHHHHHCCC'
                        #EEEEEEEEEEEEEEEEeBEEbBEEbEEE

sequence              = 'TIDQWLLKNAKEDAIAELKKAGITSDFYFNAINKAKTVEEVNALKNEILKAHA'
ss_sequence           = 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'
                         #TIDQWLLKNAKEDAIAELKKAGITSDFYFNAINKAKTVEEVNALKNEILKAHA




new_system            = '1gab_ff03ua_AMBER_extended'


text = 'source leaprc.ff03ua\n'

#---------------- building the three lettle sequence -------------#
                                                                  #
three_letter_sequence = []                                        #
                                                                  #
n = 0                                                             #
for i in sequence:                                                #
	three_letter_sequence.append(three_letter[i])                 #
	                                                              #
                                                                  #
three_letter_sequence[0]  = "N"+ three_letter_sequence[0]         #
three_letter_sequence[-1] = "C"+ three_letter_sequence[-1]        #
                                                                  #
#building the string                                              #
                                                                  #
seq   = '{'                                                       #
                                                                  #
                                                                  #
for  i in three_letter_sequence:                                  #
	seq  =  seq + i + ' '                                         #
                                                                  #
seq   = seq + '}'                                                 #
#print seq                                                        #
                                                                  #
text = text + 'system  =  sequence '                              #
                                                                  #
text = text + seq                                                 #
text = text + '\n\n'                                              #
                                                                  #
#-----------------------------------------------------------------#

'''

#------------------------------------------------------------------------------------------------------------------------------------#
#                                                                                                                                    #
#  											 Building the phi and psi angles                                                         #
#                                                                                                                                    #
#------------------------------------------------------------------------------------------------------------------------------------#
                                                                                                                                     #
helix  = [-47,   -57, 180]                                                                                                           #
sheet  = [-139, -135, 180]                                                                                                           #
                                                                                                                                     #
coil   = [180,   180, 180]                                                                                                           #
                                                                                                                                     #
n = 1                                                                                                                                #
for i in ss_sequence:                                                                                                                #
	                                                                                                                                 #
	if i == 'C':    
		#coil[0]   = random.uniform(-180, 180)	
		#coil[1]   = random.uniform(-180, 180)
		text = text + 'impose system  {' +str(n)+ '} {{"N" "CA" "C" "N" ' +str(coil[0])+ '} {"C" "N" "CA" "C" '+str(coil[1])+'}}\n'  #
	if i == 'H':                                                                                                                     #
		text = text + 'impose system  {' +str(n)+ '} {{"N" "CA" "C" "N" ' +str(helix[0])+ '} {"C" "N" "CA" "C" '+str(helix[1])+'}}\n'#
	if i == 'E':                                                                                                                     #
		text = text + 'impose system  {' +str(n)+ '} {{"N" "CA" "C" "N" ' +str(sheet[0])+ '} {"C" "N" "CA" "C" '+str(sheet[1])+'}}\n'#
	n = n + 1                                                                                                                        #
                                                                                                                                     #
#------------------------------------------------------------------------------------------------------------------------------------#
'''


text = text + '\n\n'
text = text + 'saveamberparm system '+new_system+'.prmtop '+new_system+'.inpcrd \n'
text = text + 'savepdb system '+new_system+'.pdb'

text = text + '\nquit'
 
print text                                                                                                                           #
arq_out = open("leaprc", "w")		# creates a leaprc file
arq_out.writelines(text) 			# writes the XYZ file header

os.system('tleap leaprc')


'''
# NAMD Config file - autogenerated by NAMDgui plugin
# Author: Jan Saam,  saam@charite.de

# input
amber                   yes
ambercoor               /home/labio/Dropbox/cref_bachega/thiago_1fme.inpcrd
parmfile                thiago_1fme.prmtop

# output
set output              /home/labio/Dropbox/cref_bachega/thiago_1fme
outputname              $output
dcdfile                 ${output}.dcd
xstFile                 ${output}.xst
dcdfreq                 50
xstFreq                 50

binaryoutput            no
binaryrestart           no
outputEnergies          100
restartfreq             1000

fixedAtoms              off
gbis                    yes

# Basic dynamics
exclude                 scaled1-4
1-4scaling              1
COMmotion               no
dielectric              1.0

# Simulation space partitioning
switching               on
switchdist              9
cutoff                  10
pairlistdist            12

# Multiple timestepping
firsttimestep           0
timestep                1
stepspercycle           20
nonbondedFreq           2
fullElectFrequency      4

# Temperature control

set temperature         298
temperature             $temperature;  # initial temperature

# Scripting

minimize            1000
reinitvels          $temperature
run                 100000
'''











