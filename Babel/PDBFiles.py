from Bio.PDB import PDBParser

def load_PDB_to_system(filename):
    parser = PDBParser()
    structure = parser.get_structure('X', filename)

    for model in structure:
        c = 1
        for chain in model:
            molecule = Molecule(id=1,  name="protein")
            n = 1
            r = 1
            for pdb_residue in chain:
                residue = Residue(id=r,  name=pdb_residue.resname)
                for pdb_atom in pdb_residue:

                    atom = Atom(id=n,
                                name=pdb_atom.name,
                                pos=pdb_atom.coord)
                    n += 1

                    residue.atoms.append(atom)
                molecule.residues.append(residue)
                r += 1
    return molecule

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
                chainID = 'X'
                resi = str(residue_i.id)
                x = float(atom_i.pos[0])
                y = float(atom_i.pos[1])
                z = float(atom_i.pos[2])
                occ = 1.00
                tpF = 1.00
                segID = 'P1'
                element = atom_i.name[0]

                #text +=  "ATOM     1  " + atom +   " " +resn+ "  {:4d}    {:8.3f}{:8.3f}{:8.3f}  1.00  0.00          Na+\n".format(resi, float(k), float(i), float(j))

                #ATOM, idx, Aname, " ",resn, ' ', chainID,resi,   x,     y,     z,    occ,   tpF,        segID,element," "
                text += "{:<6s}{:5d} {:<4s}{:1s}{:<3s}{:1s}{:4s}{:<2s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}      {:<4s}{:2s}  {:2s}\n".format(ATOM,
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

