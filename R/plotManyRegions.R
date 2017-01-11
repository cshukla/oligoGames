
plotAnnoTrack <- function(gr, annoTrack, elemNames=TRUE) {
  ## check may need to be modified
  if(!all(sapply(annoTrack, function(xx) is(xx, "GRanges"))))
    stop("all elements in 'annoTrack' needs to be 'GRanges'")
  plot(start(gr), 1, type = "n", xaxt = "n", yaxt = "n", bty = "n",
       ylim = c(0, length(annoTrack) + 0.5), xlim = c(start(gr), end(gr)), xlab = "", ylab = "")
  lapply(seq(along = annoTrack), function(ii) {
    jj <- length(annoTrack) + 1- ii
    ir <- subsetByOverlaps(annoTrack[[ii]], gr)
    start(ir) <- pmax(start(ir), start(gr))
    end(ir) <- pmin(end(ir), end(gr))
    
    if(length(ir) > 0){  
      rect(start(ir)-0.5, jj - 0.15, end(ir), jj + 0.15, #col = "grey60",
           #density=30, 
           col=alpha("black", 0.1), border = "black")
      if(elemNames & ii==2){
        lastPos <- rep(NA, length(ir)-1)
        for(k in 1:length(ir)){
          irk <- ir[k,]
          rwidth <- end(gr)-start(gr)
          gwidth <- min(end(irk), end(gr)) - 
            max(start(irk), start(gr))
          textPos <- max(start(irk), start(gr)) + gwidth/2
          if (sum(!is.na(lastPos))>0){
            separation <- (textPos-lastPos[k-1])/rwidth
            if(abs(separation) <= 0.25 & k < 3){
              textPos <- textPos + sign(separation+1)*0.20*rwidth
            }else{ 
              separation <- min(abs((textPos-lastPos)/rwidth), na.rm=TRUE)
              if(abs(separation) <= 0.25){
                textPos <- lastPos[k-1] + 0.20*rwidth
              }		   
            }
          }
          lastPos[k] <- textPos		  
          text(textPos, jj-0.65, 
               labels=mcols(irk)$name, cex=1.25)
        }
      }
    }
    mtext(names(annoTrack)[ii], side = 2, at = jj, las = 1, line = 1)
  })
}


.bsPlotTitle <- function(gr, extend, main, mainWithWidth, fdr=NULL) {
  if(is.data.frame(gr))
    gr <- data.frame2GRanges(gr)
  if(length(gr) > 1) {
    warning("plotTitle: gr has more than one element")
    gr <- gr[1]
  }
  plotChr <- as.character(seqnames(gr))
  plotRange <- c(start(gr), end(gr))
  regionCoord <- sprintf("%s: %s - %s", plotChr, 
                         format(plotRange[1], big.mark = ",", scientific = FALSE),
                         format(plotRange[2], big.mark = ",", scientific = FALSE))
  if(mainWithWidth) {
    regionWidth <- sprintf("width = %s, extended = %s", 
                           format(width(gr), big.mark = ",", scientific = FALSE),
                           format(extend, big.mark = ",", scientific = FALSE))
    if(!is.null(fdr)){
      regionFDR <- sprintf("FDR: %s", format(fdr, big.mark=",", scientific=FALSE))
      regionCoord <- sprintf("%s (%s), %s", regionCoord, regionWidth, regionFDR)
    }else{
      regionCoord <- sprintf("%s (%s)", regionCoord, regionWidth)
    }    
  }else{
    if(!is.null(fdr)){
      regionFDR <- sprintf("qval: %s", format(fdr, big.mark=",", scientific=FALSE))
      regionCoord <- sprintf("%s, %s", regionCoord, regionFDR)
    }	
  }
  if(main != "") {
    main <- sprintf("%s\n%s", main, regionCoord)
  } else {
    main <- regionCoord
  }
  main
}

.bsPlotLines <- function(x, y, col, lty, lwd, plotRange) {
  if(sum(!is.na(y)) <= 1)
    return(NULL)
  xx <- seq(from = plotRange[1], to = plotRange[2], length.out = 500)
  yy <- approxfun(x, y)(xx)
  lines(xx, yy, col = col, lty = lty, lwd = lwd)
}


.bsPlotPoints <- function(x, y, z, col, pointsMinCov, maxCov, regionWidth, weightPoints, meth, index=sampIdx) {
  ptSize <- z/maxCov + 0.2
  if(!weightPoints){
    ptSize <- 0.5
  }
  points(x[z>pointsMinCov], y[z>pointsMinCov], col = col, pch = 16, cex = ptSize)
  deg <- 2
  spn <- 0.4
  if (sum(z>=pointsMinCov) < 30) {
    spn <- 0.70
    
    if (sum(z>=pointsMinCov) < 20){
      deg <- 1
    }
  }
  if (meth){
    deg <- 1
    spn <- 0.4
  }
  loess_fit <- loess(y[z>=pointsMinCov] ~ x[z>=pointsMinCov], 
                    span=spn, degree=deg)
  lines(x[z>=pointsMinCov], predict(loess_fit), col = col, lty=index)
}


# plotting functions hacked to allow non-smoothed regions
.plotSmoothData <- function(BSseq, region, extend, addRegions, col, lty, lwd, regionCol,
                            addTicks, addPoints, pointsMinCov, highlightMain, smoothed,
                            weightPoints, meth) {
  gr <- bsseq:::.bsGetGr(BSseq, region, extend)
  BSseq <- subsetByOverlaps(BSseq, gr)
  
  ## Extract basic information
  sampleNames <- sampleNames(BSseq)
  names(sampleNames) <- sampleNames
  positions <- start(BSseq)
  if(smoothed){ smoothPs <- getMeth(BSseq, type = "smooth") }
  rawPs <- getMeth(BSseq, type = "raw")
  coverage <- getCoverage(BSseq)
  
  if (!meth){
    rawPs = getCoverage(BSseq, type = "M")
    coverage = getCoverage(BSseq, type = "Cov") - rawPs
    rawPs <- log2( (rawPs + 0.5) / (coverage + 0.5))
    coverage = getCoverage(BSseq, type = "Cov")
  }
  
  ## get col, lwd, lty
  colEtc <- bsseq:::.bsGetCol(object = BSseq, col = col, lty = lty, lwd = lwd)
  
  ## The actual plotting
  ylabel <- "Methylation"
  ylimit <- c(0,1)
  
  if (!meth){
    ylabel <- "log2 ratio (to total)"
    ylimit <- range(rawPs)
  }
  plot(positions[1], 0.5, type = "n", xaxt = "n", #yaxt = "n",
       ylim = ylimit, xlim = c(start(gr), end(gr)), xlab = "", ylab = ylabel)
  
  #axis(side = 2, at = c(0.2, 0.5, 0.8))
  if(addTicks)
    rug(positions)
  
  if(is.list(addRegions) & !is.data.frame(addRegions)){
    if(length(addRegions)>2){
      stop("Only two sets of regions can be highlighted")
    }
    if(length(regionCol)==1){
      regionCol <- c(regionCol, alpha("blue", 0.2))
      #message("Second addRegion highlight color missing.  Using transparent blue.")
    }
    bsseq:::.bsHighlightRegions(regions = addRegions[[1]], gr = gr, ylim = ylimit,
                                regionCol = regionCol[1], highlightMain = highlightMain)
    bsseq:::.bsHighlightRegions(regions = addRegions[[2]], gr = gr, ylim = ylimit,
                                regionCol = regionCol[2], highlightMain = highlightMain)
    
  }else{
    bsseq:::.bsHighlightRegions(regions = addRegions, gr = gr, ylim = ylimit,
                                regionCol = regionCol, highlightMain = highlightMain)		
  }
  
  # if(addPoints) {
  #     sapply(1:ncol(BSseq), function(sampIdx) {
  #         abline(v = positions[rawPs[, sampIdx] > 0.1], col = "grey80", lty = 1)
  #     })
  # } # This adds vertical grey lines so we can see where points are plotted
  
  if (smoothed){
    if(!addPoints){
      sapply(1:ncol(BSseq), function(sampIdx) {
        bsseq:::.bsPlotLines(positions, smoothPs[, sampIdx], col = colEtc$col[sampIdx],
                             lty = colEtc$lty[sampIdx], lwd = colEtc$lwd[sampIdx],
                             plotRange = c(start(gr), end(gr)))
      })
    }
  }else{
    if(!addPoints){
      sapply(1:ncol(BSseq), function(sampIdx) {
        bsseq:::.bsPlotLines(positions, rawPs[, sampIdx], col = colEtc$col[sampIdx],
                             lty = colEtc$lty[sampIdx], lwd = colEtc$lwd[sampIdx],
                             plotRange = c(start(gr), end(gr)))
      })}
  }
  if(addPoints) {
    sapply(1:ncol(BSseq), function(sampIdx) {
      .bsPlotPoints(positions, rawPs[, sampIdx], coverage[, sampIdx],
                    col = colEtc$col[sampIdx], pointsMinCov = pointsMinCov,
                    maxCov=max(coverage), regionWidth=end(gr)-start(gr),
                    weightPoints=weightPoints, meth=meth, index=sampIdx)
    })
  }
}

plotRegion <- function(BSseq, region = NULL, extend = 0, main = "", addRegions = NULL, 
                       annoTrack = NULL, col = NULL, lty = NULL, lwd = NULL,
                       BSseqTstat = NULL, stat = "tstat.corrected", stat.col = "black",
                       stat.lwd = 1, stat.lty = 1, stat.ylim = c(-8,8),
                       mainWithWidth = TRUE,
                       regionCol = alpha("red", 0.2), addTicks = TRUE, addPoints = FALSE,
                       pointsMinCov = 5, highlightMain = FALSE, smoothed=TRUE, fdr=NULL,
                       newPars=TRUE, weightPoints=TRUE, meth=TRUE) {
  if (newPars){
    if(is.null(BSseqTstat))
      layout(matrix(1:2, ncol = 1), heights = c(2,1))
    else
      layout(matrix(1:3, ncol = 1), heights = c(2,2,1))
  }
  
  .plotSmoothData(BSseq = BSseq, region = region, extend = extend, addRegions = addRegions,
                  col = col, lty = lty, lwd = lwd, regionCol = regionCol,
                  addTicks = addTicks, addPoints = addPoints,
                  pointsMinCov = pointsMinCov, highlightMain = highlightMain, smoothed=smoothed,
                  weightPoints=weightPoints, meth=meth)
  gr <- bsseq:::.bsGetGr(BSseq, region, extend)
  
  if(!is.null(BSseqTstat)) {
    BSseqTstat <- subsetByOverlaps(BSseqTstat, gr)
    plot(start(gr), 0.5, type = "n", xaxt = "n", yaxt = "n",
         ylim = stat.ylim, xlim = c(start(gr), end(gr)), xlab = "", ylab = "t-stat")
    axis(side = 2, at = c(-5,0,5))
    abline(h = 0, col = "grey60")
    mapply(function(stat, col, lty, lwd) {
      bsseq:::.bsPlotLines(start(BSseqTstat), BSseqTstat@stats[, stat],
                           lty = lty, plotRange = c(start(gr), end(gr)), col = col, lwd = lwd)
    }, stat = stat, col = stat.col, lty = stat.lty, lwd = stat.lwd)
  }
  
  if(!is.null(main)) {
    if (!is.null(fdr)){
      fdr <- round(fdr, 4)
      main <- .bsPlotTitle(gr = region, extend = extend, main = main,
                           mainWithWidth = mainWithWidth, fdr=fdr)
    }else{
      main <- .bsPlotTitle(gr = region, extend = extend, main = main,
                           mainWithWidth = mainWithWidth)
    }
    mtext(side = 3, text = main, outer = FALSE, cex = 1, line=0 )
  }
  
  if(!is.null(annoTrack))
    plotAnnoTrack(gr, annoTrack)
  
  abline(h=0, col="grey", lty=2, lwd=2)
  
  # return(invisible(NULL))
}

# modify plot many regions plot to take in a vector of fdr values (1 per region in regions 
# argument) to be displayed in the plot title.

# set addPoints = TRUE to plot individual points sized by coverage and one smooth (loess)
# line per sample instead of a uniform-sized verbatim line going through each observation

#			plotManyRegions(bs, regions=dmrs.ols.only, smoothed=FALSE, 
# 				extend=(dmrs.ols.only$end - dmrs.ols.only$start + 1)/2, 
# 				addRegions=dmrs.ols.only, annoTrack = annot, FDR=dmrs.ols.only$qval,
# 				addPoints=TRUE, pointsMinCov=1)

#			plotManyRegions(bs, regions=dmrs.ols.only, smoothed=FALSE, 
# 				extend=(dmrs.ols.only$end - dmrs.ols.only$start + 1)/2, 
# 				addRegions=dmrs.ols.only, annoTrack = annot, FDR=dmrs.ols.only$qval,
# 				addPoints=TRUE, pointsMinCov=1)

plotManyRegions <- function(BSseq, regions = NULL, extend = 0, main = "", addRegions = NULL,
                            annoTrack = NULL, col = NULL, lty = NULL, lwd = NULL, 
                            BSseqTstat = NULL, stat = "tstat.corrected", stat.col = "black",
                            stat.lwd = 1, stat.lty = 1, stat.ylim = c(-8,8),
                            mainWithWidth = TRUE, regionCol = alpha("red", 0.1), 
                            addTicks = TRUE, addPoints = FALSE, pointsMinCov = 5, 
                            highlightMain = FALSE, verbose = TRUE, smoothed=TRUE, 
                            FDR=NULL, newPars=TRUE, weightPoints=TRUE, meth=TRUE) {
  cat("[plotManyRegions] preprocessing ...")
  if(!is.null(regions)) {
    if(is(regions, "data.frame"))
      gr <- data.frame2GRanges(regions, keepColumns = FALSE)
    else
      gr <- regions
    if(!is(gr, "GRanges"))
      stop("'regions' needs to be either a 'data.frame' (with a single row) or a 'GRanges' (with a single element)")
  } else {
    gr <- granges(BSseq)
  }
  gr <- resize(gr, width = 2*extend + width(gr), fix = "center")
  BSseq <- subsetByOverlaps(BSseq, gr)
  if(!is.null(BSseqTstat))
    BSseqTstat <- subsetByOverlaps(BSseqTstat, gr)
  
  if(length(start(BSseq)) == 0)
    stop("No overlap between BSseq data and regions")
  if(!is.null(main) && length(main) != length(gr))
    main <- rep(main, length = length(gr))
  cat("done\n")
  if (newPars){
    opar <- par(mar = c(0,4.1,0,0), oma = c(1,0,3,1))
    on.exit(par(opar))
  }
  if (length(extend)==1){
    extend <- rep(extend, length(gr))
  }
  for(ii in seq(along = gr)) {
    if(verbose & ii%%100==0) cat(sprintf("[plotManyRegions]   plotting region %d (out of %d)\n", ii, nrow(regions)))
    if(!is.null(FDR)){
      fdr <- round(FDR[ii], 4) 
    }else{
      fdr <- NULL
    }
    
    
    plotRegion(BSseq = BSseq, region = regions[ii,], extend = extend[ii],
               col = col, lty = lty, lwd = lwd, main = main[ii], BSseqTstat = BSseqTstat,
               stat = stat, stat.col = stat.col, stat.lwd = stat.lwd,
               stat.lty = stat.lty, stat.ylim = stat.ylim,
               addRegions = addRegions, regionCol = regionCol, mainWithWidth = mainWithWidth,
               annoTrack = annoTrack, addTicks = addTicks, addPoints = addPoints,
               pointsMinCov = pointsMinCov, highlightMain = highlightMain,
               smoothed=smoothed, fdr=fdr, newPars=newPars, weightPoints=weightPoints, meth=meth)
  }
}

alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}
