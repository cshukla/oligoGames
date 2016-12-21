# 'DMR-finder' code adapted to find regions in lncRNAs that have differential
# oligo reporter signal in nucleus versus cytosol signal

# here maxGap parameter is meaningless since the gap between neighboring sites
# is always one (sequencing)
# can consider each lncRNA as a separate 'chromosome'
# want to increase minimum number of sites per differential region (~20?)

##Loading the Necessary Libraries

#' A function to
#'
#' This function performs
#' @param dat
#' @param minNum
#' @param bpSpan
#' @param minInSpan
#' @param verbose
#' @import GenomicRanges
#' @import bumphunter
#' @import bsseq
#' @import reshape2
#' @import dplyr
#' @import tidyr
#' @import ggplot2
#' @import splines
#' @import doParallel
#' @keywords inference
#' @return list of fitted clusters smoothened by locfit

locfitByCluster0 <- function(dat=NULL,
                             minNum = 7, bpSpan = 1000, minInSpan = 0,
                             verbose = TRUE) {
  ## Depends on limma
  ## bpSpan is in basepairs
  ## assumed x are ordered
  ## if y is vector change to matrix
  if(is.null(dat))
    stop("Data input needed for smoother")

  y = matrix(dat$y, ncol=1)
  x = dat$x
  weights = matrix(dat$weights, ncol=1)
  cluster = dat$cluster


  if(is.null((dat$y)))
    stop("y (rawBeta) is missing")
  if(is.null(dat$x))
    stop("x (pos) is missing")
  if(is.null(dat$cluster))
    stop("cluster is missing")
  if(is.null(dat$weights))
    weights <- matrix(1, nrow = nrow(y), ncol = ncol(y))
  rm(dat)

  Indexes <- split(seq(along = cluster), cluster)
  clusterL <- sapply(Indexes, length)
  smoothed <- rep(TRUE, nrow(y))

  for(i in seq(along=Indexes)) {
    if(verbose) if(i %% 1e4 == 0) cat(".")
    Index <- Indexes[[i]]
    if(clusterL[i] >= minNum & sum(rowSums(is.na(y[Index,,drop=FALSE]))==0) >= minNum)  {
      nn <- minInSpan / length(Index)
      for(j in 1:ncol(y)) {
        sdata <- data.frame(pos = x[Index],
                            y = y[Index, j],
                            weights = weights[Index,j])
        fit <- locfit(y ~ lp(pos, nn = nn, h = bpSpan), data = sdata,
                      weights = weights, family = "gaussian", maxk = 10000)
        pp <- preplot(fit, where = "data", band = "local",
                      newdata = data.frame(pos = x[Index]))
        y[Index,j] <- pp$trans(pp$fit)
      }
    } else {
      y[Index,] <- NA
      smoothed[Index] <- FALSE
    }
  }
  return(list(fitted=y, smoothed=smoothed, smoother="locfit"))
}

#' A function to
#'
#' This function performs
#' @param x
#' @param f
#' @param cutoff Defaults to quantile(abs(x), 0.99)
#' @param assumeSorted Defaults to FALSE
#' @param verbose Defaults to FALSE
#' @import
#' @keywords inference
#' @return

getSegments <- function(x, f = NULL, cutoff=quantile(abs(x), 0.99),
                        assumeSorted = FALSE, verbose=FALSE){
  if(is.null(f))
    f <- rep(1L, length(x))
  stopifnot(length(x) == length(f))
  stopifnot(length(cutoff) <= 2)
  if(is.character(f))
    f <- as.factor(f)
  if(is.numeric(f))
    f <- as.integer(f)
  stopifnot(is.factor(f) || is.integer(f))
  if(length(cutoff) == 1)
    cutoff <- c(-cutoff, cutoff)
  cutoff <- sort(cutoff)

  reordered <- FALSE
  if(!assumeSorted && is.unsorted(f)) {
    od <- order(f)
    x <- x[od]
    f <- f[od]
    reordered <- TRUE
  }

  if(verbose) message("getSegments: segmenting")
  Indexes <- split(seq(along=x), f)
  direction <- as.integer(greaterOrEqual(x, cutoff[2]))
  direction[x <= cutoff[1]] <- -1L

  ## We now need to convert f into cid
  if(verbose) message("getSegments: splitting")
  segments <- cumsum(c(1, diff(direction) != 0)) +
    cumsum(c(1, diff(f) != 0))
  names(segments) <- NULL

  res <- list(upIndex = split(which(direction>0), segments[direction>0]),
              dnIndex = split(which(direction<0), segments[direction<0]),
              zeroIndex = split(which(direction==0), segments[direction==0]))
  names(res[[1]]) <- NULL
  names(res[[2]]) <- NULL
  names(res[[3]]) <- NULL

  if(reordered) {
    res <- lapply(res, function(sp) lapply(sp, function(xx) od[xx]))
  }
  res
}

#' A function to
#'
#' This function performs
#' @param chr
#' @param pos
#' @param maxGap Defaults to 300
#' @param assumeSorted Defaults to FALSE
#' @import
#' @keywords inference
#' @return list of segments

clusterMaker <- function(chr, pos, assumeSorted = FALSE, maxGap=300){
  nonaIndex <- which(!is.na(chr) & !is.na(pos))
  Indexes <- split(nonaIndex, chr[nonaIndex])
  clusterIDs <- rep(NA, length(chr))
  LAST <- 0
  for(i in seq(along = Indexes)){
    Index <- Indexes[[i]]
    x <- pos[Index]
    if(!assumeSorted){
      Index <- Index[order(x)]
      x <- pos[Index]
    }
    y <- as.numeric(diff(x) > maxGap)
    z <- cumsum(c(1, y))
    clusterIDs[Index] <- z + LAST
    LAST <- max(z) + LAST
  }
  clusterIDs
}

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
#' @param maxGap Defaults to 300
#' @param cutoff Defaults to quantile(abs(x), 0.99)
#' @param assumeSorted Defaults to FALSE
#' @param meth.mat Defaults to meth.mat
#' @param unmeth.mat Defaults to unmeth.mat
#' @param verbose Defaults to TRUE
#' @param design Defaults to design
#' @param workers Defaults to workers
#' @param betabin
#' @import
#' @keywords inference
#' @return

# change function so that it simultaneously calculates statistic(s) while
# finding regions...
regionFinder <- function(x, chr, pos, cluster=NULL, y=x,
                         ind=seq(along=x),order=TRUE, minNumRegion=5,
                         maxGap=300, cutoff=quantile(abs(x), 0.99),
                         assumeSorted = FALSE, meth.mat=meth.mat,
                         unmeth.mat = unmeth.mat, verbose = TRUE,
                         design=design, workers=workers, betabin){
  if(any(is.na(x[ind]))){
    warning("NAs found and removed. ind changed.")
    ind <- intersect(which(!is.na(x)),ind)
  }
  if(is.null(cluster))
    cluster <- clusterMaker(chr, pos, maxGap=maxGap, assumeSorted = assumeSorted)
  Indexes <- getSegments(x = x[ind], f = cluster[ind], cutoff = cutoff,
                         assumeSorted = assumeSorted, verbose = verbose)
  clusterN <- table(cluster)[as.character(cluster)]

  # only keep up and down indices
  Indexes <- Indexes[1:2]

  # only keep clusters that have at least minNumRegion elements
  for(i in 1:2){
    lns <- sapply(Indexes[[i]], length)
    Indexes[[i]] <- Indexes[[i]][lns >= minNumRegion]
  }


  stat.lncRNA.spline <- function(ix){
    X = as.vector(sapply(design[,2], function(x) rep(x, nrow(meth.mat[ix,]))))
    L = as.vector(rep(pos[ix], nrow(design)))
    Y = log2( (meth.mat[ix,] + 1) / (unmeth.mat[ix,] + 1) )
    Y <- melt(Y)

    Y$X <- X
    Y$L <- L
    int.knots <- ceiling((max(pos[ix]) - min(pos[ix])) / 100) + 1
    Y$wts <- melt(log2(meth.mat[ix, ] + unmeth.mat[ix, ] + 1))$value
    return(summary(lm(value~ X + ns(L, df=int.knots), data=Y, weights=wts))$coef[2,3])
  }


  # perform weighted least squares version to weight each
  # observation by the total number of reads mapping to each oligo
  stat.lncRNA.loci <- function(ix){
    X = as.vector(sapply(design[,2], function(x) rep(x, nrow(meth.mat[ix,]))))
    L = as.vector(rep(pos[ix], nrow(design)))
    Y = log2( (meth.mat[ix,] + 1) / (unmeth.mat[ix,] + 1) )
    Y <- melt(Y)

    Y$X <- X
    Y$L <- L
    Y$wts <- melt(log2(meth.mat[ix, ] + unmeth.mat[ix, ] + 1))$value
    return(summary(lm(value~ X + factor(L), data=Y, weights=wts))$coef[2,3])
  }


  # stat.betabin function
  stat.betabin.fn <- function(ix){
    df <- data.frame(
      g.fac=factor(as.vector(sapply(design[,2], function(x)
        rep(x, nrow(meth.mat[ix,]))))),
      meth=melt(meth.mat[ix,])$value,
      unmeth=melt(unmeth.mat[ix,])$value
    )
    fit1 <- tryCatch({betabin(formula = cbind(meth, unmeth) ~ g.fac, ~ 1,
                              data=df, link = "logit", method="BFGS")},
                     error=function(cond) {return(NULL)},
                     warning=function(cond) {return(NULL)})

    fit2 <- tryCatch({betabin(formula = cbind(meth, unmeth) ~ 1, ~ 1,
                              data=df, link = "logit", method="BFGS")},
                     error=function(cond) {return(NULL)},
                     warning=function(cond) {return(NULL)})
    if (length(fit1)==0 | length(fit2)==0){
      fit1 <- tryCatch({betabin(formula = cbind(meth, unmeth) ~ g.fac, ~ 1,
                                data=df, link = "logit", method="Nelder-Mead")},
                       error=function(cond) {return(NULL)},
                       warning=function(cond) {return(NULL)})

      fit2 <- tryCatch({betabin(formula = cbind(meth, unmeth) ~ 1, ~ 1,
                                data=df, link = "logit", method="Nelder-Mead")},
                       error=function(cond) {return(NULL)},
                       warning=function(cond) {return(NULL)})
    }

    if (!(length(fit1)==0 | length(fit2)==0)){
      return(as.numeric(anova(fit2, fit1)@anova.table[2,9]))
    }else{
      return(NA)
    }
  }

  stat.wls.spline.paired.fn <- function(ix){
    X = as.vector(sapply(design[,2], function(x) rep(x, nrow(meth.mat[ix,]))))
    Y = meth.mat[ix,] / (meth.mat[ix,] + unmeth.mat[ix,])
    Y = Y[,design[,2]==1] - Y[,design[,2]==0]
    rownames(Y) = pos[ix]
    Y <- melt(Y)
    wt.mat <- meth.mat[ix,design[,2]==1] + unmeth.mat[ix,design[,2]==1] +
      meth.mat[ix,design[,2]==0] + unmeth.mat[ix,design[,2]==0]
    Y$wts <- melt(wt.mat)$value
    int.knots <- ceiling((max(pos[ix]) - min(pos[ix])) / 500) + 1
    return(summary(lm(value~ ns(Var1, df=int.knots), weights=wts, data=Y))$coef[1,3])
  }

  cov.upper <- quantile(meth.mat + unmeth.mat, 0.9)
  t1 <- proc.time()
  res <- vector("list",2)
  for(i in 1:2){
    res[[i]]<-
      data.frame(chr=sapply(Indexes[[i]], function(Index) chr[ind[Index[1]]]),
                 start=sapply(Indexes[[i]], function(Index) min(pos[ind[Index]])),
                 end=sapply(Indexes[[i]], function(Index) max(pos[ind[Index]])),
                 value=sapply(Indexes[[i]], function(Index) mean(y[ind[Index]])),
                 area=sapply(Indexes[[i]], function(Index)abs(sum(y[ind[Index]]))),
                 cluster=sapply(Indexes[[i]], function(Index)cluster[ind[Index]][1]),
                 indexStart=sapply(Indexes[[i]], function(Index) min(ind[Index])),
                 indexEnd = sapply(Indexes[[i]], function(Index) max(ind[Index])),
                 L = sapply(Indexes[[i]], length), stringsAsFactors=FALSE)

    if(paired){
      res[[i]]$stat.ols = unlist(mclapply(Indexes[[i]], function(Index)
        stat.wls.spline.paired.fn(ind[Index]),
        mc.cores=workers))
    }else{
      res[[i]]$stat.ols = unlist(mclapply(Indexes[[i]], function(Index)
        stat.lncRNA.spline(ind[Index]),
        mc.cores=workers))
    }
    if (betabin){
      #res[[i]]$stat.betabin = unlist(mclapply(Indexes[[i]], function(Index)
      #                 			stat.betabin.fn(ind[Index]),
      #                 			mc.cores=workers))

      # sub in weighted least squares estimator for beta binomial
      # since beta binomial is too time and computationally-expensive to
      # carry out and wls is a faster alternative that is able to consider
      # differences in coverage which result in differences in our confidence
      # in methylation proportion estimates at individual loci
      res[[i]]$stat.betabin = unlist(mclapply(Indexes[[i]], function(Index)
        stat.lncRNA.spline(ind[Index]),
        mc.cores=workers))

    }else{
      res[[i]]$stat.betabin <- NA
    }

    #res[[i]]$L <- res[[i]]$indexEnd - res[[i]]$indexStart+1
    #res[[i]]$clusterL <- sapply(Indexes[[i]], function(Index) clusterN[ind[Index][1]])

  }
  t2 <- proc.time()
  t2 - t1
  names(res) <- c("up","dn")
  res <- rbind(res$up,res$dn)
  if(order & nrow(res)>0) res <- res[order(-res$area),]

  return(res)
}

#' A function to
#'
#' This function performs
#' @param x
#' @param y
#' @import
#' @keywords inference
#' @return

greaterOrEqual <- function(x,y) {
  precision <- sqrt(.Machine$double.eps)
  (x >= y) | (abs(x-y) <= precision)
}


#' A function to
#'
#' This function performs
#' @param y
#' @param x Defaults to NULL
#' @param weights Defaults to NULL
#' @param chr
#' @param maxGapSmooth
#' @param smoothFunction
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
                     smoothFunction, verbose=TRUE, workers=1, minNum = minNum,
                     minInSpan = minInSpan, bpSpan = bpSpan) {

  #' A function to
  #'
  #' This function performs
  #' @param ix
  #' @import
  #' @keywords inference
  #' @return

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
    y <- matrix(y, ncol=1) ##need to change this to be more robust
  if(!is.null(weights) && is.null(dim(weights)))
    weights <- matrix(weights, ncol=1)
  if (!getDoParRegistered())
    registerDoSEQ()

  ret.all <- NULL
  cores <- workers
  # loop through each chromosome to mitigate memory spikes
  for (chromosome in unique(chr)){

    clusterC <- clusterMaker(chr[chr==chromosome], x[chr==chromosome], maxGap = maxGapSmooth)
    while (length(unique(clusterC)) < cores){
      newGap <- maxGapSmooth*0.9
      clusterC <- clusterMaker(chr[chr==chromosome],
                               x[chr==chromosome], maxGap = newGap)
      maxGapSmooth <- newGap
    }


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
    ret <- foreach(idx = iter(IndexesChunks), .packages = "bumphunter") %dorng% {
      sm <- locfitByCluster2(idx)
      c(sm, list(idx = seq(1, length(chr))[chr==chromosome][idx]))
    }

    attributes(ret)[["rng"]] <- NULL
    ## Paste together results from different workers
    ret <- reduceIt(ret)

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

#' A function to
#'
#' This function performs
#' @param x
#' @param elem
#' @param bind Defaults to rbind
#' @import
#' @keywords inference
#' @return

reduceIt <- function(x, elem, bind=rbind) {
  if(missing(elem)) elem <- names(x[[1]])
  ret <- lapply(elem, function(el) {
    xx <- lapply(x, "[[", el)
    if (is.matrix(xx[[1]])) return(do.call(bind, xx))
    else if (is.vector(xx[[1]])) return(do.call("c", xx))
    else stop("reduce can only handle matrices or vectors")
  })
  names(ret) <- elem
  ret
}

#' A function to
#'
#' This function performs
#' @param meth.mat
#' @param unmeth.mat
#' @param design
#' @param coef
#' @import
#' @keywords inference
#' @return

getEstimatePooled = function(meth.mat, unmeth.mat, design, coef){
  rowSums(meth.mat[,design[,coef]==1]) /
    (rowSums(meth.mat[,design[,coef]==1] + unmeth.mat[,design[,coef]==1])) -
    rowSums(meth.mat[,design[,coef]==0]) /
    (rowSums(meth.mat[,design[,coef]==0] + unmeth.mat[,design[,coef]==0]))
}

#' A function to
#'
#' This function performs
#' @param mat
#' @param wt.mat
#' @param design
#' @param coef
#' @import
#' @keywords inference
#' @return

getEstimateWeighted = function (mat, wt.mat, design, coef)
{
  designw <- apply(design, 2, function(x) x * sqrt(wt.mat))
  matw <- mat * sqrt(wt.mat)

  fit <- lm(mat ~ design[,2], weights=wt.mat)
  a <- coef(summary(fit))[1,1]
  b <- coef(summary(fit))[2,1]
  sigma <- summary(fit)$sigma
  vcov.mat <- vcov(fit)/(sigma^2)
  A1 = exp(a + b)/(1 + exp(a + b))^2 - exp(a)/(1 + exp(a))^2
  A2 = exp(a + b)/(1 + exp(a + b))^2
  var.coef = A1^2*vcov.mat[1, 1] + 2*A1*A2*vcov.mat[1, 2] + A2^2*vcov.mat[2, 2]
  meth.diff = exp(a + b)/(1 + exp(a + b)) - exp(a)/(1 + exp(a))
  sd.meth.diff = sqrt(var.coef)*sigma
  if (sigma==0){
    sd.meth.diff <- 0
  }
  out <- c(meth.diff, sd.meth.diff)
  return(out)
}

#' A function to compute the coefficient estimates for regression of the
#' methylation levels on the group indicator variable, at each CpG site
#'
#' This function computes the coefficient estimates for regression of the
#' methylation levels on the group indicator variable, at each CpG site
#'
#' @param mat
#' @param design
#' @param coef
#' @param B Defaults to NULL
#' @param permutations Defaults to NULL
#' @param full Defaults to FALSE
#' @import
#' @keywords inference
#' @return

# Function to compute the coefficient estimates for regression of the
# methylation levels on the group indicator variable, at each CpG site
getEstimate = function (mat, design, coef)
{
  v <- design[, coef]
  A <- design[, -coef, drop = FALSE]
  qa <- qr(A)
  S <- diag(nrow(A)) - tcrossprod(qr.Q(qa))
  vv <- matrix(v, ncol = 1)
  
  sv <- S %*% vv
  vsv <- diag(crossprod(vv, sv))
  b <- (mat %*% crossprod(S, vv)/vsv)
  a <- mat %*% as.matrix(rep(1/nrow(design), nrow(design))) - b*mean(v) 

  x.cov = design[, 2]
  deno = sum((x.cov - mean(x.cov))^2)
  vcov.mat = matrix(c(mean(x.cov^2)/deno, rep(-mean(x.cov)/deno, 2), 1/deno), 2, 2)
  A1 = exp(a + b)/(1 + exp(a + b))^2 - exp(a)/(1 + exp(a))^2
  A2 = exp(a + b)/(1 + exp(a + b))^2
  var.coef = A1^2*vcov.mat[1, 1] + 2*A1*A2*vcov.mat[1, 2] + A2^2*vcov.mat[2, 2]
  meth.diff = exp(a + b)/(1 + exp(a + b)) - exp(a)/(1 + exp(a))

  if (!is.matrix(b))
    b <- matrix(b, ncol = 1)
  if (!is.matrix(mat))
    mat <- matrix(mat, nrow = 1)

  sy <- mat %*% S
  df.residual <- ncol(mat) - qa$rank - 1
  sigma <- matrix(sqrt(rowSums((sy - tcrossprod(b,
                                                sv))^2)/df.residual), ncol = 1)
  sd.meth.diff = sqrt(var.coef)*sigma
  out <- list(meth.diff = meth.diff, sd.meth.diff = sd.meth.diff)
 
  return(out)
}

#' A function for the WGBS version of the '.getEstimate' function in Bumphunter
#'
#' This function is the WGBS version of the '.getEstimate' function in
#' Bumphunter for computing the estimated difference of population means for the
#' two biological groups
#'
#' @param meth.mat
#' @param unmeth.mat
#' @param design
#' @param method Defaults to LM
#' @param coef
#' @param B Defaults to NULL
#' @param permutations Defaults to NULL
#' @param full
#' @param workers Defaults to NULL
#' @import
#' @keywords inference
#' @return

##The WGBS version of the '.getEstimate' function in Bumphunter for computing
##the estimated difference of population means for the two biological groups
estim = function (meth.mat, unmeth.mat, design, coef)
{
  mat <- log2( (meth.mat + 0.5) / (unmeth.mat + 0.5) )
  
  lm.fit = getEstimate(mat = mat, design = design, coef = coef)
  meth.diff = lm.fit$meth.diff
  sd.meth.diff = lm.fit$sd.meth.diff
 
  ## Returning the list of estimated mean methylation differences and their SD s at all CpG sites
  out <- list(meth.diff = meth.diff, sd.meth.diff = sd.meth.diff)
  return(out)
}

#' A function for the WGBS version of the 'BumphunterEngine' function with 'LM'
#' and 'GLM'
#'
#' This function is the WGBS version of the 'BumphunterEngine' function with 'LM'
#' and 'GLM'
#'
#' @param meth.mat
#' @param unmeth.mat
#' @param design
#' @param method Defaults to LM
#' @param chr Defaults to NULL
#' @param pos
#' @param cluster Defaults to NULL
#' @param coef
#' @param minInSpan Defaults to 70
#' @param minNum Defaults to 70
#' @param minNumRegion Defaults to 5
#' @param cutoff Defaults to NULL
#' @param pickCutoff Defaults to FALSE
#' @param pickCutoffQ Defaults to 0.99
#' @param maxGap Defaults to 500
#' @param maxGapSmooth Defaults to 1e8
#' @param nullMethod Defaults to GLMM
#' @param smooth Defaults to FALSE
#' @param bpSpan Defaults to 1000
#' @param betabin Defaults to FALSE
#' @param smoothFunction locfitByCluster
#' @param useWeights Defaults to TRUE
#' @param B Defaults to ncol(permutations)
#' @param permutations Defaults to NULL
#' @param verbose Defaults to TRUE
#' @param workers Defaults to NULL
#' @param loci Defaults to TRUE
#' @param subject Defaults to TRUE
#' @import
#' @keywords inference
#' @return


# set input parameters
workers <- 1 # don't parallelize 
sampleSize <- 4
meth.mat = meth.mat
unmeth.mat = unmeth.mat
design = design
coef = 2
cutoff = 0.05
maxGap = 500
maxGapSmooth = 50
bpSpan = 100
minNum = 3
minNumRegion = 3
minInSpan = 4 # 10
smooth = FALSE
smoothFunction = locfitByCluster
useWeights = TRUE
verbose = TRUE

maxPerms = 34



##The WGBS version of the 'BumphunterEngine' function with 'LM' and 'GLM'
wgbs.bumphunter = function (meth.mat, unmeth.mat, design, 
                            chr = NULL, pos, coef = 2, minInSpan=70, minNum=70, minNumRegion=5,
                            cutoff = NULL, maxGap = 500, maxGapSmooth=1e8,
                            smooth = FALSE, bpSpan=1000,
                            smoothFunction = locfitByCluster, useWeights = TRUE, 
                            verbose = TRUE, workers=NULL, loci=TRUE, subject=TRUE, ...)
{
  if (!is.matrix(meth.mat))
    stop("'meth.mat' and 'unmeth.mat' must be a matrices.")
  if (ncol(meth.mat) != nrow(design))
    stop("Total number of columns in 'meth.mat' and 'unmeth.mat' must  match number of rows of 'design'")
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
  
  cluster <- clusterMaker(chr, pos, maxGap = maxGap)
  
  if (verbose)
    message("Computing coefficients.")

  tmp = estim(meth.mat = meth.mat, unmeth.mat = unmeth.mat, design = design,
              coef = coef)
  rawBeta <- tmp$meth.diff
  
  if (useWeights & smooth) {
    sd.raw = tmp$sd.meth.diff
    weights <- sd.raw 
    weights[sd.raw < 10^-5] = mean(sd.raw[sd.raw > 10^-5]) 
    rm(sd.raw)
  }
  rm(tmp)
  gc()

  if (smooth) {
    if (verbose)
      message("Smoothing coefficients")

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
                      meth.mat = meth.mat, unmeth.mat = unmeth.mat, 
                      verbose = TRUE,
                      design = design, workers=workers,  betabin=betabin)
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
#' @param meth.mat
#' @param unmeth.mat
#' @param design
#' @param method Defaults to LM
#' @param chr Defaults to NULL
#' @param pos
#' @param cluster Defaults to NULL
#' @param coef
#' @param minInSpan Defaults to 10
#' @param minNum Defaults to 70
#' @param minNumRegion Defaults to 5
#' @param cutoff Defaults to NULL
#' @param pickCutoff Defaults to FALSE
#' @param pickCutoffQ Defaults to 0.99
#' @param maxGap Defaults to 500
#' @param maxGapSmooth Defaults to 1e8
#' @param nullMethod Defaults to GLMM
#' @param smooth Defaults to FALSE
#' @param bpSpan Defaults to 1000
#' @param betabin Defaults to FALSE
#' @param smoothFunction locfitByCluster
#' @param useWeights Defaults to TRUE
#' @param B Defaults to ncol(permutations)
#' @param permutations Defaults to NULL
#' @param verbose Defaults to TRUE
#' @param workers Defaults to NULL
#' @param loci Defaults to TRUE
#' @param subject Defaults to TRUE
#' @param sampleSize Defaults to ncol(meth.mat)/2
#' @import
#' @keywords inference
#' @return
#' @export

DRfinder <- function(lncRNA, minInSpan=10,
                      design, chr = NULL, pos, 
                      coef = 2, minNum=70, bpSpan=1000, maxGapSmooth=1e8, minNumRegion=5, 
                      cutoff = NULL, maxGap = 500,
                      smooth = FALSE, smoothFunction = locfitByCluster,
                      useWeights = TRUE,  verbose = TRUE,
                      workers = NULL, loci=TRUE, subject=TRUE, sampleSize=(ncol(meth.mat)/2),
                      maxPerms=50){

  # create BSSeq object out of the lncRNA counts (for ease in plotting fns 
  #  already set up)
  # M = the nuclear or cytosol signal
  # Cov = the total signal plus the nuclear or cytosol signal 
  #  (have to constrain 0<=M<=Cov)
  
  bs <- BSseq(chr = rownames(lncRNA), pos = lncRNA[,1],
              M = lncRNA[,2:9], Cov = lncRNA[,c(10:13, 10:13)] + lncRNA[,2:9], 
              sampleNames = colnames(lncRNA)[2:9])
  pData(bs)$Type <- c(rep("Cytosol", 4), rep("Nuclear", 4))
  
  design = model.matrix(~c(rep("Cytosol", 4), rep("Nuclear", 4)))
  
  nlocs = nrow(lncRNA) 
  meth.mat = lncRNA[,2:9]
  unmeth.mat = lncRNA[,c(10:13, 10:13)]
  
  chr = rownames(lncRNA)
  pos = lncRNA[,1]
  colnames(meth.mat) <- colnames(unmeth.mat) <-  sampleNames(bs)
  
  # get observed stats

  res = wgbs.bumphunter(meth.mat = meth.mat, unmeth.mat = unmeth.mat, minInSpan=minInSpan,
                        design = design, chr = chr, pos = pos, 
                        coef = coef, minNum=minNum, maxGapSmooth=maxGapSmooth, minNumRegion=minNumRegion,
                        cutoff = cutoff, maxGap = maxGap, 
                        smooth = smooth, smoothFunction = smoothFunction, 
                        useWeights = useWeights,  verbose = verbose,
                        workers = workers, loci=loci, subject=subject) 


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
      
      res.flip.p = wgbs.bumphunter(meth.mat = meth.mat, unmeth.mat = unmeth.mat, 
                                   minInSpan=minInSpan,
                                   design = designr, chr = chr, pos = pos, 
                                   coef = coef, minNum=minNum, maxGapSmooth=maxGapSmooth, minNumRegion=minNumRegion,
                                   cutoff = cutoff,  maxGap = maxGap, 
                                   smooth = smooth, smoothFunction = smoothFunction, 
                                   useWeights = useWeights,  verbose = verbose,
                                   workers = workers, loci=loci, subject=subject) 
      res.flip.p$permNum <- j 
      
      res.flip <- rbind(res.flip, res.flip.p)
      print(paste0(j, " out of ", ncol(perms), " permutations completed"))
    }	
  }else{
    stop("Error: Currently only balanced designs supported")
  }
  rm(res.flip.p)
  
  perm.ordered <- c(sort(abs(res$stat.ols), method="quick"), Inf)
  pval <- rep(NA, nrow(res))
  pval[!is.na(res$stat.ols)] <- (1 + 
                          sapply(abs(res$stat.ols[!is.na(res$stat.ols)]), 
                                function(x) length(perm.ordered) - 
                                min(which(x <= perm.ordered)))) /
                           (1 + sum(!is.na(res$stat.ols)))							
  res$pval <- pval
  res$qval <- p.adjust(pval, method="BH")
  return(res)
}
