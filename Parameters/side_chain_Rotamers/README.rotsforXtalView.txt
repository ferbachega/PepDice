README.rotsforXtalView.txt 
Modification Date: 21 November 01 (wba) 
Creation Date: 13 June 00 (scl) 
Web Ref: http://kinemage.biochem.duke.edu/databases/rotamer.html


INTRODUCTION

All-atom contact analysis shows that all published rotamer
libraries contain serious van der Waals overlaps [1].  This
should not occur as rotamers, being the more common
conformations, should have the lower energy states.  Using a
select database of 240 high resolution, low-clash score [2], low
Rcryst structures and then filtering it by B-factor and clash
score, we have composed a rotamer library, consisting of 153
conformers, which we think is more faithful to the rotamer
concept and will improve accuracy of new structures. The library
is available as an O database, and this document discussed the
use of that database.


USE INSTRUCTIONS

Download the tar file from our web or ftp site.

The rotamer libraries are in a tar'ed Zipped file named 
rotsForXView.tar.Z.  The tar file can be downloaded via your
browser from either the web site - rotamer section or by
maneuvering through the FTP site hierarchy: 
pub/datasets/rotamers. Or, if you prefer, use "anonymous" login
to the FTP site.

1. Extract the files from the tar file.

These commands will work: 
   uncompress rotsForXView.tar.Z 
   tar -xvf 
   rotsForXView.tar cd rotsForXView

2. Copy the database to the XtalView data directory.

First create a backup copy of the current rotamer libraries,
dict.noh.pdb and dict.allh.pdb.  Then, copy/move the two new
rotamer dictionaries, into the XtalView data directory (typically
/usr/local/XtalView/data/  ).

You can now select either of these dictionaries from the "Models"
pane of Xfit.


COMMENTS

The rotamer dictionary  was made from a database of 240
structures at 1.7A resolution or better. Care was taken to remove
side chains which had uncertain positions (eg high B-factors) or
systematically misfit conformations. The citation below [1]
describes our methodology, and the advantages of this library
over others previously published.

The approach taken in preparing this file differs from that of
the standard rotamer file.

Firstly, although the most common rotamer comes first, subsequent
rotamers are arranged so that neighbours have similar
conformations. This allows you to narrow down your choice to a
few rotamers which will be next to each other as you turn the
dial.

Secondly, it contains many more rotamers. We examined all
rotamers and distributions to avoid (we hope) artifacts.  Thus
there are some rare-but-real conformations.  This means that
non-rotameric conformations should be viewed with even more
suspicion than before.

Thirdly, it contains some non-rotamers! These are sample points
in allowed, well populated, relaxed conformations.  They are not,
however, local energy minima and thus, strictly, not rotamers.
For example, Glu chi3 has a very flat distribution which is
populated throughout, so we have included the rotamer at chi3=0
and two sample points at chi3=60 and chi3=-60. There are sample
points only for Glu, Gln, Asp, Asn, Phe and Tyr. They all have
their frequency listed as 0%.

KNOWN BUGS: Histidine can have several protonation states.  The
allh file only allows for one (doubly protonated).  This
behaviour is the same as the original files provided with
XtalView.


REFERENCES

[1] SC Lovell, JM Word, JS Richardson and DC Richardson (2000)
"The Penultimate Rotamer Library", Proteins: Structure, Function
and Genetics, 40: 389-408

[2]	Word, et. al. (1999) "Visualizing and Quantifying
Molecular Goodness-of-fit: Small-probe Contact Dots with Explicit
Hydrogen Atoms", J. Mol. Biol. 285: 1711-1733

Please send questions and comments to us rather than Duncan
McRee.

e-mail: Bryan Arendall (arendall@duke.edu) or David C. Richardson
(dcr@kinemage.biochem.duke.edu) URL:
http://kinemage.biochem.duke.edu Biochemistry Department Duke
University Durham, NC USA 27710



