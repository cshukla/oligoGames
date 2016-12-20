#!/usr/bin/env python

from optparse import OptionParser
import sys, os, gzip

########################################################################
# mapToBarcodes.py
#
# Map sequencing reads to oligo barcodes and generate a table of 
# counts in case and control condition for all oligos. 
#
# Notes:
# 	- Ensure the barcode exactly matches the read (NO MISMATCH ALLOWED)
# 	- Allow up to 2 mismatches in the region upstream the barcode
# 	- Compile counts for each transcript in a simple table with each
# 	  replicate as a column and individual transcripts as rows
########################################################################
def main():
	usage='usage:%prog [options] <fastqCase1,fastqCase2,fastqCase3,...> <fastqControl1,fastqControl2,fastqControl3,...>'
	parser = OptionParser(usage)
	parser.add_option('--labels', dest='labels', default='nuclear,total', type=str, help='Labels for the case and control (input) fastq files [Default: %default]')
	parser.add_option('--oligoMap', dest='str', default='oligoMap.tsv', help='A file mapping each transcript (oligo) to a unique barcode and variable sequence [Default: %default]')
	parser.add_option('--out', dest='outDir', default='oligoOut', type=str, help='Directory for output files [Default: %default]')
	(options, args) = parser.parse_args()

	if len(args)==2:
		fastqCase = args[0].split(',')
		fastqControl = args[0].split(',')
	else:
		parser.error('Please provide fastq files for case and control' %(usage))

	################################################
	# Read in the counts file headers and infer the
	# number of replicates and set the names of the
	# two conditions
	################################################
	print >> sys.stderr, 'Reading and setting input arguments...'
	labels = options.labels.split(',')
	oligoMap = options.oligoMap

	################################################
	# Check the presence of the output directory.
	# If absent create it
	################################################
	print >> sys.stderr, 'Output will be stored in %s...' %(options.outDir)
	if os.path.exists(options.outDir):
		outDir = options.outDir
	else:
		os.mkdir(options.outDir)
		outDir = options.outDir

	#######################################
	# Make a hash linking each barcode,
	# variable sequence and transcript ID
	#######################################
	print >> sys.stderr, 'Generating a hash mapping each barcode to a unique variable sequence and transcriptID'
	barcodeTranscript = {}
	for line in open(oligoMap, 'r'):
		a = line.strip().split('\t')
		transcriptID = a[0]
		barcode = a[1]
		varSeq = a[2]
		if barcode not in barcodeTranscript:
			barcodeTranscript[barcode] = {}
		barcodeTranscript[barcode][varSeq] = transcriptID

	################################################
	# Compute the # of times every barcode appears
	# in the set of input FASTQ files.
	################################################
	print >> sys.stderr, 'Mapping the fastq files to %s\n Output will be written to %s' %(oligoSeqs, allOligoCounts)
	print >> sys.stderr, 'This step is expected to be slow so please be patient - get a coffee or something'
	allOligoCounts = '/'.join([outDir, 'allTranscriptsCounts.tsv'])
	barcodeCounts(fastqCase, fastqControl, labels, allOligoCounts, barcodeTranscript)
	inputData = allOligoCounts

######################################################################
# barcodeCounts
#
# Input:
# 	- fastqCase: List of FASTQ files for case (condition 1)
# 	- fastq_control: List of FASTQ files for control (condition 2)
# 	- labels: A list of labels for case and control (nuclear, total)
# 	- allOligoCounts: Output file to write counts of all oligos
# 	- oligoSeqs: Sequences of each oligo in the oligo pool
# 	- oligoInfo: File with information about the oligo pool
#
# Output:
# 	- Write counts for all oligos in each case and control fastq file
# 	  to allOligoCounts file
######################################################################
def barcodeCounts(fastqCase, fastqControl, labels, allOligoCounts, barcodeTranscript):
	print >> sys.stderr, '\n'
	# Generate header for the output file
	header = ['Transcript']

	for i in range(0, len(fastqCase)):
		label = labels[0] + str(i+1)
		header.append(label)

	for i in range(0, len(fastqControl)):
		label = labels[1] + str(i+1)
		header.append(label)

	# Open the output file and write the header
	outFile = open(allOligoCounts, 'w')
	outLine = '\t'.join(header)
	print >> outFile, outLine

	barcodes = barcodeTranscript.keys()

	# Initiliaze case barcodes counts
	#print >> sys.stderr, 'Initializing hash of case and control barcode counts with empty lists'
	caseCounts = {}
	for barcode in barcodes:
		varSeq = barcodeTranscript[barcode]
		for seq in varSeq:
			if barcode not in caseCounts:
				caseCounts[barcode] = {}
			caseCounts[barcode][seq] = []

	# Initialize control barcodes counts
	controlCounts = {}
	for barcode in barcodes:
		varSeq = barcodeTranscript[barcode]
		for seq in varSeq:
			if barcode not in controlCounts:
				controlCounts[barcode] = {}
			controlCounts[barcode][seq] = []

	#######################################
	# Find counts for each barcode in all
	# the case and control fastq files.
	#######################################

	#print >> sys.stderr, 'Finding counts for each barcode in input fastq files.'
	for fastq in fastqCase:
		print >> sys.stderr, 'Processing file %s. Please wait...' %(fastq)
		fastqFile = gzip.open(fastq, 'r')
		sequences = []
		count = 3
		for line in fastqFile:
			a = line.strip()
			if count % 4 == 0:
				sequences.append(a)
			count += 1

		fastqCounts = {}

		for seq in sequences:
			barcodeRC = seq[0:10]
			varSeqRC = seq[10:]
			barcode = reverseComplement(barcodeRC)
			varSeq = reverseComplement(varSeqRC)
			if barcode not in fastqCounts:
				fastqCounts[barcode] = {}
			fastqCounts[barcode][varSeq] = 1 + fastqCounts.get(barcode, {}).get(varSeq, 0)
		
		for barcode in caseCounts:
			varSeq = caseCounts[barcode]
			readVar = fastqCounts.get(barcode, {})
			if len(varSeq) == 1:
				seqCounts = 0
				for var in readVar:
					seqCounts += readVar.get(var, 0)
				caseCounts[barcode][varSeq.keys()[0]].append(seqCounts)
			else:
				for seq in varSeq:
					seqCounts = 0
					for var in readVar:
						if hammingDistance(var, seq) <= 2:
							seqCounts += 1
					caseCounts[barcode][seq].append(seqCounts)

	for fastq in fastqControl:
		print >> sys.stderr, 'Processing file %s. Please wait...' %(fastq)
		fastqFile = gzip.open(fastq, 'r')
		sequences = []
		count = 3
		for line in fastqFile:
			a = line.strip()
			if count % 4 == 0:
				sequences.append(a)
			count += 1

		fastqCounts = {}

		for seq in sequences:
			barcodeRC = seq[0:10]
			varSeqRC = seq[10:]
			barcode = reverseComplement(barcodeRC)
			varSeq = reverseComplement(varSeqRC)
			if barcode not in fastqCounts:
				fastqCounts[barcode] = {}
			fastqCounts[barcode][varSeq] = 1 + fastqCounts.get(barcode, {}).get(varSeq, 0)
		
		for barcode in controlCounts:
			varSeq = controlCounts[barcode]
			readVar = fastqCounts.get(barcode, {})
			if len(varSeq) == 1:
				seqCounts = 0
				for var in readVar:
					seqCounts += readVar.get(var, 0)
				controlCounts[barcode][varSeq.keys()[0]].append(seqCounts)
			else:
				for seq in varSeq:
					seqCounts = 0
					for var in readVar:
						if hammingDistance(var, seq) <= 2:
							seqCounts += 1
					controlCounts[barcode][seq].append(seqCounts)

	#######################################
	# Write to output the counts of each
	# transcript in all the case and 
	# control replicates 
	#######################################

	print >> sys.stderr, 'Finished processing all the fastq files' 
	print >> sys.stderr, '\n'
	print >> sys.stderr, 'Printing counts for each transcript to %s' %(outFile)
	for barcode in barcodeTranscript:
		transcripts = barcodeTranscript[barcode]
		for seq in transcripts:
			lncrna, oligoNum = transcripts[seq]
			case = caseCounts[barcode][seq]
			control = controlCounts[barcode][seq]
			transcriptName = '_'.join([lncrna, oligoNum, barcode])
			outLine = [transcriptName] + case + control
			for i in range(0, len(outLine)):
				outLine[i] = str(outLine[i])

			outLine = '\t'.join(outLine)
			print >> outFile, outLine

	outFile.close()
	print >> sys.stderr, 'Done'

#########################################################################################
# reverseComplement
#
# Input:
# 	- seq ---> Sequence String
#
# Output:
# 	- bases ---> Reverse Complement of Sequence String
#########################################################################################
def reverseComplement(seq):
	altMap = {'ins':'0'}
	complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
	for k,v in altMap.iteritems():
		seq = seq.replace(k,v)

	bases = list(seq) 
	bases = reversed([complement.get(base,base) for base in bases])
	bases = ''.join(bases)
	for k,v in altMap.iteritems():
		bases = bases.replace(v,k)
	return bases

#########################################################################################
# hammingDistance
#
# Input:
# 	- s1: First string
# 	- s2: Second string
#
# Output:
# 	- hamming distance (# of mismatches) between the two input strings.
#########################################################################################
def hammingDistance(s1, s2):
	if len(s1) != len(s2):
		raise ValueError("Hamming Distance undefined for sequences of unequal length")
	return sum(bool(ord(ch1) - ord(ch2)) for ch1, ch2 in zip(s1, s2))