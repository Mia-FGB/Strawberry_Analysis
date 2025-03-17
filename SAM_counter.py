#!/usr/bin/python

import sys, getopt, errno
fastqFilename = "sample.fq"
samFilename = "unique_alignments.sam"
refFilename = "refs.txt"

class AlignmentCount:
    def __init__(self, ref_name):
        self.ref_name = ref_name
        self.count = 0
        self.bp = 0

try:
	opts, args = getopt.getopt(sys.argv[1:],"f:s:r:",["ifile="])
except getopt.GetoptError():
	print "Option not recognised."
	sys.exit(2)
for opt, arg in opts:
	if opt in ("-f"):
		fastqFilename = arg
	elif opt in ("-s"):
		samFilename = arg
	elif opt in ("-r"):
		refFilename = arg

readLengths = dict()	
counts = dict()	
refs = dict()

try:
	with open(refFilename, 'r') as refFile:
		for line in refFile:
			line = line.strip()
			fields = line.split()
			refName = fields[1].strip()
			seqName = fields[0].strip()
			if not refName in refs:
				refs[refName] = []
			refs[refName].append(seqName)
except (OSError, IOError) as e: 
	if getattr(e, 'errno', 0) == errno.ENOENT:
		print "Could not find file " + refFilename
	sys.exit(2)

try:
	with open(fastqFilename, 'r') as fastqFile:
		count = 0
		readName = ""
		for line in fastqFile:
			count += 1
			if count % 4 == 1:
				line = line.strip()
				fields = line.split()
				readName = fields[0][1:]
			elif count % 4 == 2:
				line = line.strip()
				readLengths[readName] = len(line)
except (OSError, IOError) as e: 
	if getattr(e, 'errno', 0) == errno.ENOENT:
		print "Could not find file " + fastqFile
	sys.exit(2)

try:
	with open(samFilename, 'r') as samFile:
		for line in samFile:
			line = line.strip()
			fields = line.split()
			readName = fields[0].strip()
			seqName = fields[2].strip()
			readLength = readLengths[readName]
			refName = ""

			for key in refs:
				if seqName in refs[key]:
					refName = key
					break

			if not refName:
				print "Error, could not find ref for sequence " + fields[2].strip()
				continue

			if refName not in counts:
				counts[refName] = AlignmentCount(refName)
			counts[refName].count += 1
			counts[refName].bp += readLength
except (OSError, IOError) as e: 
	if getattr(e, 'errno', 0) == errno.ENOENT:
		print "Could not find file " + samFilename
		sys.exit(2)

for key in counts:
	alignmentCount = counts[key]
	print alignmentCount.ref_name + "\t" + str(alignmentCount.count) + "\t" + str(alignmentCount.bp)
