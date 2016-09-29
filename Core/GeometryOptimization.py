from CRDFiles import save_CRD_to_file, load_CRD_from_file
import subprocess
import tempfile


def minimize(molecule=None,
             log=False,
             imin=1,
             maxcyc=1000,
             ncyc=100,
             cut=10,
             rgbmax=999,
             igb=1,
             ntb=0,
             ntpr=100,
             ntr=0,
             restraintmask=':1-50&@CA,N,C,O=',
             restraint_wt=50.0,
             ):

    if molecule.ff_type == 'amber':

        with tempfile.NamedTemporaryFile() as minimization_file:
            with tempfile.NamedTemporaryFile() as coordinates_file:
                with tempfile.NamedTemporaryFile() as output_file:
                    write_AMBER_minimization_input_file(
                        molecule=molecule,
                        filename=minimization_file.name,
                        imin=imin,
                        maxcyc=maxcyc,
                        ncyc=ncyc,
                        cut=cut,
                        rgbmax=rgbmax,
                        igb=igb,
                        ntb=ntb,
                        ntpr=ntpr,
                        ntr=ntr,
                        restraintmask=restraintmask,
                        restraint_wt=restraint_wt
                    )

                    save_CRD_to_file(
                        molecule=molecule, filename=coordinates_file.name)

                    subprocess.check_call([
                        'sander',
                        '-O',
                        '-i', minimization_file.name,
                        '-c', coordinates_file.name,
                        '-o', 'minimized.log',
                        '-p', molecule.top,
                        '-r', output_file.name,
                        '-ref', coordinates_file.name,
                    ])
                    load_CRD_from_file(
                        molecule=molecule, filename=output_file.name)


def write_AMBER_minimization_input_file(molecule=None,
                                        filename='minimize.in',
                                        imin=1,
                                        maxcyc=1000,
                                        ncyc=100,
                                        cut=10,
                                        rgbmax=999,
                                        igb=1,
                                        ntb=0,
                                        ntpr=100,
                                        ntr=0,
                                        restraintmask=':1-50&@CA,N,C,O=',
                                        restraint_wt=50.0):

    text = """Stage 1 - minimisation of TC5b
 &cntrl
  imin=%i, maxcyc=%i, ncyc=%i,
  cut=%i, rgbmax=%i,igb=%i, ntb=%i,
  ntpr=%i,
  ntr=%i,
  restraintmask = '%s',
  restraint_wt=%3.1f,
 /
 """ % (imin, maxcyc, ncyc, cut, rgbmax, igb, ntb, ntpr, ntr,
        restraintmask, restraint_wt)

    output_file = open(filename, "w")
    output_file.write(text)
    output_file.close()
