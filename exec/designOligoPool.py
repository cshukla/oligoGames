import gzip, sys, itertools, os, random

##############################################################
# design oligo pool
#
# Things to Do:
#   1. Add sequence names when generating oligo pool
#   2. Print out the meta data files automatically
##############################################################

def designOligoPool(regionsFile, tileSize, overlap, barcodesPerSequence, barcodesFile, univPrimers, reSeq, microSeedsFile, badNucs, numScrambles, outDir):

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
        barcodes = generateBarcodes(microSeedsFile, badNucs)

    ##############################################################
    #
    ##############################################################

    print >> sys.stderr, 'Tiling the input regions'
    variableSeqs = tileRegions(regionsFile, tileSize, overlap)

    ##############################################################
    # generate oligo pool
    ##############################################################

    print >> sys.stderr, 'There are %d variable sequences of %d bp if we use %d bp overlap' %(len(variableSeqs), tileSize, overlap)

    if len(barcodes) < float(len(variableSeqs))/(barcodesPerSequence + numScrambles):
        print 'We have more variable sequences than barcodes\nCannot generate oligo pool'
        exit(1)

    oligoPool = {}
    barCount = 0
    seqCount = 1

    for seq in variableSeqs:
        ####################################
        # For each sequence, use # barcodes
        # specified by user
        ####################################
        i = 0
        for i in range(0,barcodesPerSequence): # We want to have 10 barcodes for each variable sequence
            oligoSeq = univPrimers[0] + seq + reSeq + barcodes[barCount]
            oligoHeader = '>' + 'Sequence_' + str(seqCount) + '_' + str(i)
            oligoPool[oligoHeader] = oligoSeq
            barCount += 1
            i += 1

        ####################################
        # Scramble each sequence as
        # specified by user
        ####################################
        j = 0
        for j in range(0,numScrambles):
            scrambleSeq = list(seq)
            random.shuffle(scrambleSeq)
            scrambleOligo = univPrimers[0] + ''.join(scrambleSeq) + reSeq + barcodes[seqCount]
            scrambleHeader = '>' + 'Scramble_' + str(seqCount) + '_' + str(j)
            oligoPool[scrambleHeader] = scrambleOligo
            barCount += 1
            j += 1
        
        ####################################
        # Move on to the next sequence
        ####################################
        seqCount += 1

    ##############################################################
    # print oligo pool to output file
    ##############################################################

    if os.path.exists(outDir):
        outDir = outDir
    else:
        os.mkdir(outDir)

    outFile = '/'.join([outDir, 'oligoPool.fa'])
    outFile = open(outFile, 'w')
    for seq in oligoPool.keys():
        print >> outFile, seq
        print >> outFile, oligoPool[seq]

##############################################################
# generate barcodes
##############################################################

def generateBarcodes(microSeedsFile, badStrings):
    random.seed(21191)
    nucleotides = ['A', 'C', 'G', 'T']
    allBarcodes = itertools.product(nucleotides, repeat = 10)
    allBarcodes = [''.join(i) for i in itertools.product(nucleotides, repeat = 10)]
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
    seqs = []
    discardedSeqs = 0
    for line in open(regionsFile, 'r'):
        a = line.strip()
        if a[0] != '>':
            seqs.append(a)
    
    tiledSeqs = []
    for seq in seqs:
        start = 0
        while start <= len(seq) - tileSize + 1:
            if start + tileSize <= len(seq):
                tiledSeq = seq[start:start+tileSize]
                tiledSeqs.append(tiledSeq)
            else:
                tiledSeq = seq[-tileSize:]
                tiledSeqs.append(tiledSeq)
            start += overlap
        
    return tiledSeqs