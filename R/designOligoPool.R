#' Design oligo pool by tiling a set of sequences
#' 
#' 
#' This function helps a user design an oligo pool based on their needs. Many
#' different options provided to the user allow a high degree of customization
#' and flexibility. There are 3 main steps in this function: (1). Generate a
#' potential set of barcodes avoiding the presence of any unwanted sub strings
#' such as microRNA seeds (2) Tile each sequence in the input file using the
#' specified tile size and overlap (3) Generate final oligo pool by
#' concatenating the universal primers, potential restriction enzyme sequence
#' and a unique barcode with each tile. Each oligo is designed as
#' UPS_FP----VarSeq----RESite----Barcode----UPS_RP
#' 
#' @param regionsFile the name of the file with regions to tile for the oligo 
#'   pool. Please make sure the file is in FASTA format and has unique headers 
#'   for sequence of each region.
#' @param tileSize numeric indicating the size of each individual oligo. The 
#'   tile size is also sometimes called the length of variable sequence. Since 
#'   universal primer sites, barcodes and potentially a Restriction Enzyme site 
#'   has to be added to the tile for a complete oligo, it is usually ~50 bp less
#'   than the maximum length of the oligo you can synthesize. Default value is 
#'   110.
#' @param overlap numeric describing how many nucleotides you wish to move 
#'   between 2 successive tiles. Sometimes, this is called the step size and 
#'   this number will determine the resolution of your data downstream. Default 
#'   value is 10.
#' @param barcodesPerSequence numeric for the number of unique barcodes to use 
#'   for each variable sequence. In certain applications, the user might wish to
#'   use upto 20 barcodes for each individual tile for more accurate 
#'   measurement. Default value is 1.
#' @param barcodesFile the name of the file with a list of barcodes you wish to 
#'   use for designing the oligo pool. You don't need to give this file since 
#'   the script can generate barcodes on its own but is given to the user as an 
#'   option. If you provide the file, please note \code{microSeedsFile} and 
#'   \code{badNucs} are redundant. Default value is NULL.
#' @param univPrimers character vector of the universal primer sequences added 
#'   upstream and downstream in each oligo. Provide in the order 
#'   c('ForwardPrimer', 'ReversePrimer'). Default value is 
#'   c('ACTGGCCGCTTCACTG','AGATCGGAAGAGCGTCG'). Use non-standard primer 
#'   sequences at your own risk.
#' @param reSeq character describing the restriction enzyme site you wish to 
#'   add. To play some games, a restriction enzyme site is neccessary but this 
#'   is not strictly required. Defaults to NULL.
#' @param microSeedsFile the name of a file with microRNA seeds. These sequences
#'   are required to ensure they are not present in any barcodes. Please ensure 
#'   the file is gzipped and contains one microRNA seed sequence per line. This 
#'   is a required parameter if \code{barcodeFile} is not given. Default value 
#'   is NULL.
#' @param badNucs character vector with potential sub strings you wish to avoid 
#'   in barcodes. For example, triple N's are hard to PCR and so it is a good 
#'   idea to avoid them in barcodes. This option is highly recommended if 
#'   \code{barcodeFile} is not specified. Default value is c('AAA','TTT', 'CCC',
#'   'GGG').
#' @param numScrambles numeric describing the number of times each tile 
#'   (variable sequence) should be scrambled. Default is 0.
#' @param outDir the name of a directory you want to write the oligo pool 
#'   sequences to. The directory will be created if not present already. Default
#'   value is lncTilingGame
#' @import rPython
#' @keywords design oligo pool
#' @export

# @examples
# regionsFile <- system.file('extdata', 'lncRNAs.fa' package='oligoGames')
# outDir <- 'demoOligoGame'
# designOligoPool(regionsFile)

mapToBarcodes <- function(regionsFile, tileSize=110, overlap=10, barcodesPerSequence=1, 
                          barcodesFile=NULL, 
                          univPrimers=c('ACTGGCCGCTTCACTG','AGATCGGAAGAGCGTCG'), 
                          reSeq=NULL, microSeedsFile=NULL, 
                          badNucs=c('AAA','TTT', 'CCC', 'GGG'), numScrambles=0, 
                          outDir='lncTilingGame') {
  python.load(system.file("exec", "designOligoPool.py", package = "oligoGames"))
  python.call("designOligoPool", regionsFile, tileSize, overlap, barcodesPerSequence, barcodesFile, univPrimers, reSeq, microSeedsFile, badNucs, numScrambles, outDir)
}
