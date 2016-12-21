#' A function to
#'
#' This function performs
#' @param x
#' @param chr
#' @param pos
#' @param cluster Defaults to NULL
#' @param y Defaults to x
#' @param ind
#' @param order Defaults to TRUE
#' @param minNumRegion Defaults to 5
#' @param maxGap 
#' @param cutoff Defaults to quantile(abs(x), 0.99)
#' @param assumeSorted Defaults to FALSE
#' @param oligo.mat Defaults to oligo.mat
#' @param verbose Defaults to TRUE
#' @param design Defaults to design
#' @param workers Defaults to workers
#' @import
#' @keywords inference
#' @return

# change function so that it simultaneously calculates statistic(s) while
# finding regions...
regionFinder <- function(x, chr, pos, cluster=NULL, y=x,
                         ind=seq(along=x),order=TRUE, minNumRegion=5,
                         maxGap, cutoff=quantile(abs(x), 0.99),
                         assumeSorted = FALSE, oligo.mat=oligo.mat,
                         verbose = TRUE,
                         design=design, workers=workers){
  if(any(is.na(x[ind]))){
    warning("NAs found and removed. ind changed.")
    ind <- intersect(which(!is.na(x)),ind)
  }
  if(is.null(cluster))
    cluster <- bumphunter:::clusterMaker(chr, pos, maxGap=maxGap, 
                                         assumeSorted = assumeSorted)
  Indexes <- bumphunter:::getSegments(x = x[ind], f = cluster[ind], 
                                      cutoff = cutoff, 
                                      assumeSorted = assumeSorted, 
                                      verbose = verbose)
  clusterN <- table(cluster)[as.character(cluster)]

  # only keep up and down indices
  Indexes <- Indexes[1:2]

  # only keep clusters that have at least minNumRegion elements
  for(i in 1:2){
    lns <- sapply(Indexes[[i]], length)
    Indexes[[i]] <- Indexes[[i]][lns >= minNumRegion]
  }


  stat.OligoSignal.spline <- function(ix){
    X = as.vector(sapply(design[,2], function(x) rep(x, nrow(oligo.mat[ix,]))))
    L = as.vector(rep(pos[ix], nrow(design)))
    Y = log2(oligo.mat[ix,] + 1)
    Y <- melt(Y)

    Y$X <- X
    Y$L <- L
    int.knots <- ceiling((max(pos[ix]) - min(pos[ix])) / 100) + 1
    return(summary(lm(value~ X + ns(L, df=int.knots), data=Y))$coef[2,3])
  }

  cov.upper <- quantile(oligo.mat, 0.9)
  t1 <- proc.time()
  res <- vector("list",2)
  for(i in 1:2){
    res[[i]]<-
      data.frame(chr=sapply(Indexes[[i]], function(Index) chr[ind[Index[1]]]),
                 start=sapply(Indexes[[i]], function(Index) min(pos[ind[Index]])),
                 end=sapply(Indexes[[i]], function(Index) max(pos[ind[Index]])),
                 indexStart=sapply(Indexes[[i]], function(Index) min(ind[Index])),
                 indexEnd = sapply(Indexes[[i]], function(Index) max(ind[Index])),
                 length = sapply(Indexes[[i]], length), stringsAsFactors=FALSE)

    res[[i]]$stat = unlist(mclapply(Indexes[[i]], function(Index)
        stat.OligoSignal.spline(ind[Index]),
        mc.cores=workers))
  }
  t2 <- proc.time()
  t2 - t1
  names(res) <- c("up","dn")
  res <- rbind(res$up,res$dn)
  if(order & nrow(res)>0) res <- res[order(-res$stat),]

  return(res)
}


#' A function to
#'
#' This function performs
#' @param y
#' @param x Defaults to NULL
#' @param weights Defaults to NULL
#' @param chr
#' @param maxGapSmooth
#' @param verbose Defaults to TRUE
#' @param workers Defaults to 1
#' @param minNum
#' @param minInSpan
#' @param bpSpan
#' @import
#' @keywords inference
#' @return

# modify smoother function to behave properly with obeying the number of clusters
# specified by workers
smoother <- function(y, x=NULL, weights=NULL, chr=chr, maxGapSmooth=maxGapSmooth,
                     verbose=TRUE, workers=1, minNum = minNum,
                     minInSpan = minInSpan, bpSpan = bpSpan) {


  locfitByCluster2 <- function(ix) {

    ## if y is vector change to matrix
    yi = matrix(y[ix], ncol=1)
    xi = x[ix]
    weightsi = matrix(weights[ix], ncol=1)
    clusteri = clusterC[ix]


    if(is.null((yi)))
      stop("y (rawBeta) is missing")
    if(is.null(xi))
      stop("x (pos) is missing")
    if(is.null(clusteri))
      stop("cluster is missing")
    if(is.null(weightsi))
      weightsi <- matrix(1, nrow = nrow(yi), ncol = ncol(yi))

    Indexes <- split(seq(along = clusteri), clusteri)
    clusterL <- sapply(Indexes, length)
    smoothed <- rep(TRUE, nrow(yi))

    for(i in seq(along=Indexes)) {
      #if(verbose) if(i %% 1e4 == 0) cat(".") # print the cluster index
      Index <- Indexes[[i]]
      if(clusterL[i] >= minNum & sum(rowSums(is.na(yi[Index,,drop=FALSE]))==0) >= minNum)  {
        nn <- minInSpan / length(Index)
        for(j in 1:ncol(yi)) {
          sdata <- data.frame(posi = xi[Index],
                              yi = yi[Index, j],
                              weightsi = weightsi[Index,j])
          fit <- locfit(yi ~ lp(posi, nn = nn, h = bpSpan), data = sdata,
                        weights = weightsi, family = "gaussian", maxk = 10000)
          pp <- preplot(fit, where = "data", band = "local",
                        newdata = data.frame(posi = xi[Index]))
          yi[Index,j] <- pp$trans(pp$fit)
        }
      } else {
        yi[Index,] <- NA
        smoothed[Index] <- FALSE
      }
    }
    return(list(fitted=yi, smoothed=smoothed, smoother="locfit"))
  }

  if(is.null(dim(y)))
    y <- matrix(y, ncol=1) 
  if(!is.null(weights) && is.null(dim(weights)))
    weights <- matrix(weights, ncol=1)
  if (!getDoParRegistered())
    registerDoSEQ()

  ret.all <- NULL
  cores <- workers
  # loop through each chromosome to mitigate memory spikes
  for (chromosome in unique(chr)){

    clusterC <- bumphunter:::clusterMaker(chr[chr==chromosome], 
                                          x[chr==chromosome], 
                                          maxGap = maxGapSmooth)


    Indexes <- split(seq(along=clusterC), clusterC)
    IndexesChunks <- vector("list", length = cores)
    baseSize <- length(Indexes) %/% cores
    remain <- length(Indexes) %% cores
    done <- 0L
    for(ii in 1:cores) {
      if(remain > 0) {
        IndexesChunks[[ii]] <- done + 1:(baseSize + 1)
        remain <- remain - 1L
        done <- done + baseSize + 1L
      } else {
        IndexesChunks[[ii]] <- done + 1:baseSize
        done <- done + baseSize
      }
    }


    IndexesChunks <- lapply(IndexesChunks, function(idxes) {
      do.call(c, unname(Indexes[idxes]))
    })

    idx <- NULL ## for R CMD check
    ret <- foreach(idx = iter(IndexesChunks)) %dorng% {
      sm <- locfitByCluster2(idx)
      c(sm, list(idx = seq(1, length(chr))[chr==chromosome][idx]))
    }

    attributes(ret)[["rng"]] <- NULL
    ## Paste together results from different workers
    ret <- bumphunter:::reduceIt(ret)

    if(is.null(ret.all)){
      ret.all <- ret
    }else{
      ret.all$smoother <- c(ret.all$smoother, ret$smoother)
      ret.all$fitted <- c(ret.all$fitted, ret$fitted)
      ret.all$smoothed <- c(ret.all$smoothed, ret$smoothed)
      ret.all$idx <- c(ret.all$idx, ret$idx)
    }
    rm(ret)
  }

  ## Now fixing order issue
  revOrder <- ret.all$idx
  names(revOrder) <- seq_along(ret.all$idx)
  revOrder <- sort(revOrder)
  revOrder <- as.integer(names(revOrder))

  ret.all$smoother <- ret.all$smoother[1]
  ret.all$fitted <- ret.all$fitted[revOrder,,drop=FALSE]
  ret.all$smoothed <- ret.all$smoothed[revOrder]
  ret.all$idx <- NULL

  return(ret.all)
}

#' A function for the WGBS version of the 'BumphunterEngine' function with 'LM'
#' and 'GLM'
#'
#' This function is the WGBS version of the 'BumphunterEngine' function with 'LM'
#' and 'GLM'
#'
#' @param oligo.mat
#' @param design
#' @param chr Defaults to NULL
#' @param pos
#' @param cluster Defaults to NULL
#' @param coef
#' @param minInSpan Defaults to 10
#' @param minNum Defaults to 10
#' @param minNumRegion Defaults to 5
#' @param cutoff Defaults to NULL
#' @param maxGap Defaults to 50
#' @param maxGapSmooth Defaults to 1e8
#' @param smooth Defaults to FALSE
#' @param bpSpan Defaults to 100
#' @param verbose Defaults to TRUE
#' @param workers Defaults to NULL
#' @import
#' @keywords inference
#' @return

##The WGBS version of the 'BumphunterEngine' function with 'LM' and 'GLM'
bumphunt = function(oligo.mat, design, 
                            chr = NULL, pos, coef = 2, minInSpan=10, minNum=10, 
                            minNumRegion=5,
                            cutoff = NULL, maxGap = 50, maxGapSmooth=50,
                            smooth = FALSE, bpSpan=100,
                            verbose = TRUE, workers=NULL, ...)
{
  if (!is.matrix(oligo.mat))
    stop("'oligo.mat' must be a matrices.")
  if (ncol(oligo.mat) != nrow(design))
    stop("Total number of columns in 'oligo.mat' must  match number of rows of 'design'")
  if (!(is.null(cutoff) || length(cutoff) %in% 1:2))
    stop("'cutoff' has to be either NULL or a vector of length 1 or 2")
  if (length(cutoff) == 2)
    cutoff <- sort(cutoff)

  if (!getDoParRegistered())
    registerDoSEQ()
  registerDoParallel()
  if (is.null(workers)) { workers <- 1 }
  backend <- getDoParName()
  version <- getDoParVersion()
  subverbose <- max(as.integer(verbose) - 1L, 0)
  if (verbose) {
    if (workers == 1) {
      mes <- "Using a single core (backend: %s, version: %s)."
      message(sprintf(mes, backend, version))
    }else {
      mes <- "Parallelizing using %s workers/cores (backend: %s, version: %s)."
      message(sprintf(mes, workers, backend, version))
    }
  }
  if (is.null(chr))
    chr <- rep("Unspecified", length(pos))
  if (is.factor(chr))
    chr <- as.character(chr)
  
  cluster <- bumphunter:::clusterMaker(chr, pos, maxGap = maxGap)
  
  if (verbose)
    message("Computing coefficients.")

  est <- bumphunter:::.getEstimate(log2(oligo.mat+1), design, coef, full=TRUE) 
  rawBeta <- est$coef
  
  if (smooth) {
    sd.raw = est$sigma*est$stdev.unscaled
    weights <- sd.raw 
    weights[sd.raw < 10^-5] = mean(sd.raw[sd.raw > 10^-5]) 
    rm(sd.raw)
  }
  rm(est)
  gc()

  if (smooth) {
    if (verbose)
      message("Smoothing coefficients.")

    beta <- vector("list", 3)
    beta[[1]] <- beta[[2]] <- rep(NA, length(pos))
    for (chromosome in unique(chr)){
      beta.tmp <- smoother(y = rawBeta[chr==chromosome],
                           x = pos[chr==chromosome],
                           workers=workers, chr=chr[chr==chromosome],
                           maxGapSmooth=maxGapSmooth, weights = weights,
                           minNum = minNum, minInSpan = minInSpan, bpSpan = bpSpan,
                           verbose = verbose)
      beta[[1]][chr==chromosome] <- beta.tmp[[1]]
      beta[[2]][chr==chromosome] <- beta.tmp[[2]]
    }
    beta[[3]] <- beta.tmp[[3]]
    names(beta) <- names(beta.tmp)
    rm(beta.tmp)
    Index <- which(beta$smoothed)
    beta <- beta$fitted

  }else {
    beta <- rawBeta
    Index <- seq(along = beta)
  }
  rm(rawBeta)
  gc()
  if (verbose)
    message("Finding regions.")
  tab <- regionFinder(x = beta, chr = chr, pos = pos, cluster = cluster,
                      cutoff = cutoff, ind = Index, minNumRegion = minNumRegion,
                      oligo.mat = oligo.mat,
                      verbose = TRUE,
                      design = design, workers=workers)
  rm(beta);
  rm(chr);
  rm(pos);
  gc()
  if (nrow(tab) == 0) {
    if (verbose)
      message("No candidate regions found!")
    return(NA)
  }else {
    if (verbose)
      message(sprintf("Found %s candidate regions.", nrow(tab)))
  }

  return(table = tab)
}

#' A function for the WGBS version of the 'BumphunterEngine' function with 'LM'
#' and 'GLM'
#'
#' This function is the WGBS version of the 'BumphunterEngine' function with 'LM'
#' and 'GLM'
#'
#' @param OligoSignal
#' @param conditionLabels
#' @param design
#' @param coef
#' @param minInSpan Defaults to 10
#' @param minNum Defaults to 70
#' @param minNumRegion Defaults to 5
#' @param cutoff Defaults to NULL
#' @param nullMethod Defaults to GLMM
#' @param smooth Defaults to FALSE
#' @param bpSpan Defaults to 1000
#' @param verbose Defaults to TRUE
#' @param workers Defaults to NULL
#' @param sampleSize Defaults to (ncol(OligoSignal)-1)/2
#' @param maxPerms 
#' @import
#' @keywords inference
#' @return
#' @export

DRfinder <- function(OligoSignal, 
                     conditionLabels=c("condition1", "condition2"),
                     minInSpan=10,
                      design,
                      coef = 2, minNum=70, bpSpan=1000, 
                      minNumRegion=5, 
                      cutoff = NULL, 
                      smooth = FALSE,
                      verbose = TRUE,
                      workers = NULL, 
                     sampleSize=(ncol(OligoSignal)-1)/2,
                      maxPerms=50){

  cond <- c(rep(conditionLabels[1], sampleSize), 
            rep(conditionLabels[2], sampleSize))
  design = model.matrix(~cond)
  
  nlocs = nrow(OligoSignal) 
  oligo.mat = OligoSignal[,-1]
  
  chr = rownames(OligoSignal)
  pos = OligoSignal[,1]
  rm(OligoSignal)
  
  maxGap <- maxGapSmooth <-max(diff(pos))+1
  
  # get observed stats

  res = bumphunt(oligo.mat, 
                        minInSpan=minInSpan,
                        design = design, chr = chr, pos = pos, 
                        coef = coef, minNum=minNum, maxGapSmooth=maxGapSmooth, 
                        minNumRegion=minNumRegion,
                        cutoff = cutoff, maxGap = maxGap, 
                        smooth = smooth, 
                        verbose = verbose,
                        workers = workers) 


  if (nrow(design)%%2==0){
    perms <- combn(seq(1, nrow(design)), sampleSize)
    perms <- perms[, 2:(ncol(perms)/2)]
    res.flip <- NULL
    
    if (maxPerms < ncol(perms)){
      # subset on 'balanced perms'
      if (sampleSize > 3 & sampleSize < 6){ 
        sg <- apply(perms, 2, function(x) sum(x > sampleSize))
        perms <- perms[, sg < (sampleSize-1) & sg >= 2 ]
        maxPerms <- min(maxPerms, ncol(perms))
      }else if(sampleSize >= 6){
        sg <- apply(perms, 2, function(x) sum(x > sampleSize))
        perms <- perms[, sg >= floor(sampleSize/2) & sg <= ceiling(sampleSize/2) ]
      }
      perms <- perms[,sort(sample(1:ncol(perms), maxPerms, replace=FALSE))]
    }
    
    # Now rerun on flipped designs and concatenate results
    for (j in 1:ncol(perms)){
      reorder <- perms[,j]
      designr <- design
      designr[,2] <- 0
      designr[reorder,2] <- 1
      
      res.flip.p = bumphunt(oligo.mat, 
                                   minInSpan=minInSpan,
                                   design = designr, chr = chr, pos = pos, 
                                   coef = coef, minNum=minNum, 
                                   maxGapSmooth=maxGapSmooth, 
                                   minNumRegion=minNumRegion,
                                   cutoff = cutoff,  maxGap = maxGap, 
                                   smooth = smooth,
                                   verbose = verbose,
                                   workers = workers) 
      res.flip.p$permNum <- j 
      
      res.flip <- rbind(res.flip, res.flip.p)
      print(paste0(j, " out of ", ncol(perms), " permutations completed"))
    }	
  }else{
    stop("Error: Currently only balanced designs supported")
  }
  rm(res.flip.p)
  
  perm.ordered <- c(sort(abs(res$stat), method="quick"), Inf)
  pval <- rep(NA, nrow(res))
  pval[!is.na(res$stat)] <- (1 + 
                          sapply(abs(res$stat[!is.na(res$stat)]), 
                                function(x) length(perm.ordered) - 
                                min(which(x <= perm.ordered)))) /
                           (1 + sum(!is.na(res$stat)))							
  res$pval <- pval
  res$qval <- p.adjust(pval, method="BH")
  return(res)
}
