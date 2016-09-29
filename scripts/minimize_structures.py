#!/usr/bin/env python

import argparse
import glob

# Can't use multiprocessing, since pepdice uses fixed input filenames
# and a second process overwrites input files from the first
# from multiprocessing import Pool
import os

import pepdice.Examples.PDBs.database_builder as pepdice


parser = argparse.ArgumentParser(
    description='Minimize pdb structures.')
parser.add_argument(
    'input_structure',
    help='Reference structure to be minimized')
parser.add_argument(
    '--output_dir', help='Output directory for the decoys', default='results')
args = parser.parse_args()


def minimize_structure(pdb_file, output_dir):
    input_filename = os.path.basename(pdb_file)[:-4]  # Remove extension
    base_filename = os.path.join(output_dir, input_filename)

    if not os.path.exists(base_filename + '_A_AMBER_minimized.pdb'):
        pepdice.pdb_extract_chain(
            pdbin=pdb_file,
            pdbout=base_filename + '_A.pdb',
            chain='A',
        )
        pepdice.get_sequence_from_pdb(
            pdbin=base_filename + '_A.pdb',
            seq_out=base_filename + '_A.seq',
        )
        pepdice.build_AMBER_system_from_PDB(
            pdbin=base_filename + '_A.pdb',
            basename=base_filename + '_A_AMBER',
            force_field='ff03ua.labio',
            overwrite=True,
        )
        system = pepdice.Molecule()
        system.load_PDB_to_system(filename=base_filename + '_A_AMBER.pdb')
        system.import_AMBER_parameters(
            top=base_filename + '_A_AMBER.top',
            torsions=os.path.join(
                pepdice.PEPDICE_PARAMETER, 'amber/AMBER_rotamers.dat'),
        )
        pepdice.minimize(
            molecule=system,
            imin=1,
            maxcyc=1000,
            ncyc=100,
            cut=10,
            rgbmax=999,
            igb=1,
            ntb=0,
            ntpr=100,
            ntr=0,
        )
        pepdice.save_PDB_to_file(
            system, base_filename + '_A_AMBER_minimized.pdb')


if __name__ == '__main__':
    # Create output dir if it doesnt exist
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    if os.path.isfile(args.input_structure):
        minimize_structure(args.input_structure, args.output_dir)
    else:
        input_dir = args.input_structure
        pdb_files = glob.glob(os.path.join(input_dir, '*.pdb'))
        pdb_files += glob.glob(os.path.join(input_dir, '*.ent'))

        print 'Minimizing {} PDB files'.format(len(pdb_files))
        for pdb_file in pdb_files:
            minimize_structure(pdb_file, args.output_dir)
