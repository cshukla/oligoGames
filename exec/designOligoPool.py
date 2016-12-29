import gzip, sys, itertools, os, random

##############################################################
# design oligo pool
#
# Things to Do:
#   1. Add sequence names when generating oligo pool
#   2. Print out the meta data files automatically
##############################################################

def designOligoPool(regionsFile, tileSize, overlap, barcodesPerSequence, barcodesFile, univPrimers, reSeq, microSeedsFile, badNucs, numScrambles, barcodeLen, outDir):
	##############################################################
	# get all barcodes
	##############################################################

	if os.path.exists(barcodesFile):
		barcodes = []
		for line in open(barcodesFile):
			a = line.strip()
			barcodes.append(a)
	else:
		print >> sys.stderr, 'Generating Barcodes. Please Wait...'
		barcodes = generateBarcodes(str(microSeedsFile), badNucs, barcodeLen)

	##############################################################
	# ensure that output directory exists and open files to 
	# print the oligo pool, oligo map and meta data
	##############################################################

	if os.path.exists(outDir):
		outDir = outDir
	else:
		os.mkdir(outDir)

	oligoPool = '/'.join([outDir, 'oligoPool.fa'])
	oligoPool = open(oligoPool, 'w')

	oligoMap = '/'.join([outDir, 'oligoMap.tsv'])
	oligoMap = open(oligoMap, 'w')
	print >> oligoMap, '\t'.join(['seqName','barcode', 'variableSeq'])

	##############################################################
	# tile the input regions 
	##############################################################

	print >> sys.stderr, 'Tiling the input regions'
	variableSeqs = tileRegions(regionsFile, tileSize, overlap)

	##############################################################
	# generate oligo pool
	##############################################################

	print >> sys.stderr, 'There are %d variable sequences of %d bp if we use %d bp overlap' %(len(variableSeqs), tileSize, overlap)
	totalOligos = (len(variableSeqs)*barcodesPerSequence) + (numScrambles*len(variableSeqs))
	print >> sys.stderr, 'We would require %d oligos to tile each region with %d barcodes and scramble each sequence %d times' %(totalOligos, barcodesPerSequence, numScrambles)

	if len(barcodes) < float(len(variableSeqs))/(barcodesPerSequence + numScrambles):
		print 'We have more variable sequences than barcodes\nCannot generate oligo pool'
		exit(1)

	barCount = 0
	seqCount = 1
	barcodeIndexes = {}
	numOfOligos = {}

	for seqName in sorted(variableSeqs.keys()):
		startBarcodeIndex = barCount
		regionName = '_'.join(seqName.split('_')[:-1])
		
		if regionName not in barcodeIndexes:
			barcodeIndexes[regionName] = []
		barcodeIndexes[regionName].append(startBarcodeIndex)

		if regionName not in numOfOligos:
			numOfOligos[regionName] = 0
		numOfOligos[regionName] += barcodesPerSequence

		seq = variableSeqs[seqName]

		####################################
		# For each sequence, use # barcodes
		# specified by user
		####################################
		i = 0
		for i in range(0,barcodesPerSequence): # We want to have 10 barcodes for each variable sequence
			oligoSeq = univPrimers[0] + seq + reSeq + barcodes[barCount]
			oligoHeader = '_'.join([seqName, str(i+1)])
			print >> oligoPool, oligoHeader
			print >> oligoPool, oligoSeq
			oligoMapOut = '\t'.join([oligoHeader[1:], barcodes[barCount], seq])
			print >> oligoMap, oligoMapOut
			barCount += 1
			i += 1

		####################################
		# Move on to the next sequence
		####################################
		seqCount += 1

	##############################################################
	# print meta file
	##############################################################
	oligoMeta = '/'.join([outDir, 'oligoMeta.tsv'])
	oligoMeta = open(oligoMeta, 'w')
	print >> oligoMeta, '\t'.join(['name', 'numOfOligos', 'window', 'oligoLen', 'startBarcodeIndex', 'seqLen'])

	seqLens = {}
	barCount = 0
	for line in open(regionsFile):
		a = line.strip()
		if a[0]=='>':
			regionName = a
			seqLens[regionName] = 0
		else:
			seqLens[regionName] += len(a)

	for regionName in sorted(seqLens.keys()):
		name = regionName[1:]
		numOligos = numOfOligos[regionName]
		window = overlap
		oligoLen = tileSize
		startIndex = min(barcodeIndexes[regionName])
		seqLen = seqLens[regionName]
		metaLine = '\t'.join([name, str(numOligos), str(window), str(oligoLen), str(startIndex), str(seqLen)])
		print >> oligoMeta, metaLine

'''
	####################################
	# Scramble each sequence as
	# specified by user
	####################################
	for seqName in variableSeqs:
		seq = variableSeqs[seqName]

		j = 0
		for j in range(0,numScrambles):
			scrambleSeq = list(seq)
			random.shuffle(scrambleSeq)
			scrambleOligo = univPrimers[0] + ''.join(scrambleSeq) + reSeq + barcodes[barCount]
			scrambleHeader = '_'.join(['>Scramble' + seqName[1:], str(j+1)])
			print >> oligoPool, scrambleHeader
			print >> oligoPool, scrambleOligo
			oligoMapOut = '\t'.join([scrambleHeader[1:], barcodes[barCount], ''.join(scrambleSeq)])
			print >> oligoMap, oligoMapOut
			barCount += 1
			j += 1
'''

##############################################################
# generate barcodes
##############################################################
def generateBarcodes(microSeedsFile, badStrings, barcodeLen):
	random.seed(21191)
	nucleotides = ['A', 'C', 'G', 'T']
	allBarcodes = itertools.product(nucleotides, repeat = barcodeLen)
	allBarcodes = [''.join(i) for i in itertools.product(nucleotides, repeat = barcodeLen)]
	microSeeds = []
	for line in open(microSeedsFile):
		microSeeds.append(line.strip())

	badStrings += microSeeds

	for i in badStrings:
		filteredBarcodes = [barcode for barcode in allBarcodes if any(bad in barcode for bad in badStrings)]
	
	filteredBarcodes = set(filteredBarcodes)
	finalBarcodes = [barcode for barcode in allBarcodes if barcode not in filteredBarcodes]
	
	return finalBarcodes

##############################################################
# generate variable sequences
##############################################################

def tileRegions(regionsFile, tileSize, overlap):
	seqs = {}
	discardedSeqs = 0
	for line in open(regionsFile, 'r'):
		a = line.strip()
		if a[0] == '>':
			seqName = a
			seqs[seqName] = ''
		else:
			seqs[seqName] += a
	
	tiledSeqs = {}
	for seqName in seqs.keys():
		seqIndex = 0
		seq = seqs[seqName]
		start = 0
		while start <= len(seq) - tileSize + 1:
			if start + tileSize <= len(seq):
				tiledSeq = seq[start:start+tileSize]
				tileName = '_'.join([seqName, str(seqIndex)])
				tiledSeqs[tileName] = tiledSeq
				seqIndex += 1
			else:
				tiledSeq = seq[-tileSize:]
				tileName = '_'.join([seqName, str(seqIndex)])
				tiledSeqs[tileName] = tiledSeq
				seqIndex += 1
			start += overlap
		
	return tiledSeqs
