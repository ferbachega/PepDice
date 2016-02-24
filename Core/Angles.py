#!/usr/bin/env python

"""phipsi: compute the phi and psi angles in a protein.
Supported input format is currently only PDB.
Supported output format is residue-angle list or phi-psi-count histogram."""

import logging
import math

def main():
	import sys, getopt, os.path

	# Set up default values
	try:
		logging.basicConfig(format=os.path.basename(sys.argv[0])
					+ ": %(levelname)s: %(message)s")
	except TypeError:
		logging.basicConfig()
	try:
		opts, args = getopt.getopt(sys.argv[1:], "dv")
	except getopt.GetoptError, e:
		logging.error(e)
		raise SystemExit, 1
	surfaceType = "SES"

	# Process command line flags
	for opt, val in opts:
		if opt == "-d":
			logging.getLogger().setLevel(logging.DEBUG)
		elif opt == "-v":
			logging.getLogger().setLevel(logging.INFO)

	# Identify input and output streams
	Stdin = "<standard input>"
	Stdout = "<standard output>"
	if len(args) == 2:
		inFilename = args[0]
		outFilename = args[1]
	elif len(args) == 1:
		inFilename = args[0]
		outFilename = "<standard output>"
	elif len(args) == 0:
		inFilename = "<standard input>"
		outFilename = "<standard output>"
	else:
		logging.error("Usage: %s [-dv] [ input_PDB [ output ]]"
				% sys.argv[0])
		raise SystemExit, 1

	# Read in data
	logging.info("Reading from %s" % inFilename)
	if inFilename is Stdin:
		inFile = sys.stdin
	else:
		try:
			inFile = open(inFilename, "r")
		except IOError, e:
			logging.error(e)
			raise SystemExit, 1
	try:
		p = read(inFile, inFilename)
	except IOError, msg:
		logging.error("%s: %s", inFilename, msg)
		raise SystemExit, 1
	inFile.close()

	# Convert atoms into protein and compute phi-psi values
	computePhiPsi(p)
	for c in p:
		logging.debug("%d amino acids in chain", len(c))

	# Write out data
	logging.info("Writing to %s" % outFilename)
	if outFilename is Stdout:
		outFile = sys.stdout
	else:
		try:
			outFile = open(outFilename, "w")
		except IOError, e:
			logging.error(e)
			raise SystemExit, 1
	try:
		write(outFile, outFilename, p)
	except IOError, msg:
		logging.error("%s: %s", outFilename, msg)
		raise SystemExit, 1
	outFile.close()
	raise SystemExit, 0

# Store protein backbone data and compute phi and psi angles.
#
# Protein p is a list of chains.
# Chain c is a list of amino acids.
# Each amino acid is a 8-list of [seq, chain, type, N, CA, C, phi, psi].
# N, CA, C is each a coordinate (a cartesian module Point).

def computePhiPsi(p):
	"Compute the phi and psi angles of all amino acids in all chains."
	for c in p:
		_computeChainPhiPsi(c)

def _computeChainPhiPsi(c):
	"Compute the phi and psi angles of all amino acids"
	for i in range(1, len(c)):
		c[i][6] = _computePhi(c[i - 1], c[i])
	for i in range(0, len(c) - 1):
		c[i][7] = _computePsi(c[i], c[i + 1])

def _computePhi(prev, this):
	"Compute the phi angle of this amino acid."
	return dihedral(prev[5], this[3], this[4], this[5])

def _computePsi(this, next):
	"Compute the psi angle of this amino acid."
	return dihedral(this[3], this[4], this[5], next[3])

def readPDB(f, fname):
	"Read in PDB format input file and return a Protein instance."
	p = []
	c = []
	rSeq = None
	rChain = None
	rType = None
	bb = {}		# backbone map
	saved = False	# has this amino acid been added already?
	for line in f:
		if line[:3] == "TER":
			# Chain break
			if c:
				p.append(c)
				c = []
			continue
		if line[:4] != "ATOM":
			continue
		atomName = line[12:16].strip()
		if atomName not in [ "N", "CA", "C" ]:
			# Ignore non-backbone atoms
			logging.debug("skipping atom %s" % atomName)
			continue
		resType = line[17:20]
		resChain = line[21]
		resSeq = int(line[22:26])
		if rSeq != resSeq or rChain != resChain or rType != resType:
			rSeq = resSeq
			rType = resType
			rChain = resChain
			bb = {}
			saved = False
		x = float(line[30:38])
		y = float(line[38:46])
		z = float(line[46:54])
		bb[atomName] = Point(x, y, z)
		logging.debug("saved atom %s, len=%d", atomName, len(bb))
		if len(bb) == 3 and not saved:
			# We saw an amino acid since all three
			# keys (N, CA and C) are in our dictionary
			aa = [rSeq, rChain, rType, bb["N"], bb["CA"], bb["C"],
				None, None]
			c.append(aa)
			saved = True
	if c:
		p.append(c)
	return p

#
# Input reader registry.
# "readers" is a map from file extensions to an input reader.
# When called to read a file (via module function "read"), we
# check to see if we recognize the file type by looking for
# an input reader for the file extension (eg ".pdb").  If so,
# we call the function and return its value; otherwise, we
# raise an exception and let the upper layers handle it.
#
readers = {}

def registerReader(ext, func):
	"Register a reader for a files with given extension."
	readers[ext] = func

registerReader(".pdb", readPDB)
registerReader(".PDB", readPDB)
registerReader(".ent", readPDB)
registerReader(".ENT", readPDB)

def read(f, fname):
	"Convert data from input file into a Protein instance."
	import os.path
	if fname[0] == '<':
		ext = ".pdb"
	else:
		root, ext = os.path.splitext(fname)
	try:
		func = readers[ext]
	except KeyError:
		raise IOError, "unsupported file extension: %s" % ext
	else:
		return func(f, fname)

# Output phi-psi data in one of several formats.
# Currently supported formats are:
# 	residue-phi-psi tab-separated-values
# 	phi-psi-count tab-separated-value

def writeRPP(f, fname, p, **kw):
	"Write out phi and psi angle for each residue as tab-separated values."
	for c in p:
		_writeRPPChain(f, fname, c)

def _writeRPPChain(f, fname, c):
	"Write out phi and psi for a single amino acid chain."
	for aa in c:
		seq = aa[0]
		chain = aa[1]
		type = aa[2] 
		phi = aa[6]
		psi = aa[7]
		values = [ str(seq), chain, type,
				_angleStr(phi), _angleStr(psi) ]
		print >> f, "\t".join(values)

def _angleStr(a):
	"Return string representation of an angle in degrees."
	if a is None:
		return "-"
	else:
		return "%.1f" % math.degrees(a)

def writePPC(f, fname, p, **kw):
	"Write out a Ramachadran histogram of phi-psi bins."
	# First we compute the counts
	numBins = int(kw.get("numBins", 18))	# number of bins
	histogram = {}				# key=(phi,psi) value=count
	for c in p:
		_addToHistogram(c, histogram, numBins)

	# Then we get the range of values
	minCount = 0
	maxCount = 0
	for count in histogram.itervalues():
		minCount = min(count, minCount)
		maxCount = max(count, maxCount)

	# Then we output the plot
	zero = '.'
	symbols = ":+*@"
	r = float(maxCount - minCount) / len(symbols)
	limits = [ int(round((i + 1) * r)) + minCount
			for i in range(len(symbols)) ]
	print >> f, "'%c' = 0" % zero
	lo = 0
	for i in range(len(symbols)):
		print >> f, "'%c' = (%d..%d]" % (symbols[i], lo, limits[i])
		lo = limits[i]
	for psiIndex in range(numBins - 1, -1, -1):
		line = []
		for phiIndex in range(0, numBins):
			try:
				count = histogram[(phiIndex, psiIndex)]
				for i in range(len(symbols)):
					if count <= limits[i]:
						line.append(symbols[i])
						break
				else:
					raise ValueError, \
						"count %d out of range" % count
			except KeyError:
				line.append(zero)
		print >> f, ''.join(line)

def _addToHistogram(c, histogram, numBins):
	"Add angles from chain into histogram counts."
	origin = -math.pi
	binSize = math.pi * 2 / numBins
	for aa in c:
		phi = aa[6]
		psi = aa[7]
		if phi is None or psi is None:
			continue
		phiIndex = _bin(phi, numBins, origin, binSize)
		psiIndex = _bin(psi, numBins, origin, binSize)
		key = (phiIndex, psiIndex)
		try:
			histogram[key] += 1
		except KeyError:
			histogram[key] = 1

def _bin(a, numBins, origin, binSize):
	"Return the bin index for angle a."
	n = int(round((a - origin) / binSize))
	if n < 0:
		return 0
	elif n >= numBins:
		return numBins - 1
	else:
		return n

#
# Output writer registry.
# "writers" is a map from file extensions to an output writer.
# When called to write a file (via module function "write"), we
# check to see if we recognize the file type by looking for
# an output writer for the file extension (eg ".rpp").  If so,
# we call the function and return its value; otherwise, we
# raise an exception and let the upper layers handle it.
#
writers = {}

def registerWriter(ext, func):
	"Register a writer for a files with given extension."
	writers[ext] = func

registerWriter(".rpp", writeRPP)
registerWriter(".RPP", writeRPP)
registerWriter(".ram", writePPC)
registerWriter(".RAM", writePPC)

def write(f, fname, p, **kw):
	"Write data from a Protein instance into a file."
	import os.path
	if fname[0] == '<':
		ext = ".ram"
	else:
		root, ext = os.path.splitext(fname)
	try:
		func = writers[ext]
	except KeyError:
		raise IOError, "unsupported file extension: %s" % ext
	else:
		return func(f, fname, p, **kw)

# Support Cartesian points, vectors and operators (dot products,
# cross products, etc).

def Point(x, y, z):
	return (x, y, z)

def Vector(x, y, z):
	return (x, y, z)

def length(v):
	"Return length of a vector."
	sum = 0.0
	for c in v:
		sum += c * c
	return math.sqrt(sum)

def subtract(u, v):
	"Return difference between two vectors."
	x = u[0] - v[0]
	y = u[1] - v[1]
	z = u[2] - v[2]
	return Vector(x, y, z)

def dot(u, v):
	"Return dot product of two vectors."
	sum = 0.0
	for cu, cv in zip(u, v):
		sum += cu * cv
	return sum

def cross(u, v):
	"Return the cross product of two vectors."
	x = u[1] * v[2] - u[2] * v[1]
	y = u[2] * v[0] - u[0] * v[2]
	z = u[0] * v[1] - u[1] * v[0]
	return Vector(x, y, z)

def angle(v0, v1):
	"Return angle [0..pi] between two vectors."
	cosa = dot(v0, v1) / length(v0) / length(v1)
	return math.acos(cosa)

def dihedral(p0, p1, p2, p3):
	"Return angle [0..2*pi] formed by vertices p0-p1-p2-p3."
	v01 = subtract(p0, p1)
	v32 = subtract(p3, p2)
	v12 = subtract(p1, p2)
	v0 = cross(v12, v01)
	v3 = cross(v12, v32)
	# The cross product vectors are both normal to the axis
	# vector v12, so the angle between them is the dihedral
	# angle that we are looking for.  However, since "angle"
	# only returns values between 0 and pi, we need to make
	# sure we get the right sign relative to the rotation axis
	a = angle(v0, v3)
	if dot(cross(v0, v3), v12) > 0:
		a = -a
	return a

if __name__ == "__main__":
	main()
