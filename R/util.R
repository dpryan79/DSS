###############################################
## DelayedMatrix stats functions
##
## Many operations on a DelayedMatrix object scale very poorly as its size increases. 
## This can be combatted by computing the operations on chunks of the DelayedMatrix. 
## The following functions facilitate that.
##
## A DelayedMatrix is processed in blocks of blockSize and then in parallel according to BPPARAM.
###############################################
block_mean <- function(m, na.rm=FALSE, blockSize=10000, BPPARAM=SerialParam()) {
    # mean() is VERY slow on a large delayed matrix, but quick on subset
    grid = RegularArrayGrid(refdim=dim(m), spacings=c(blockSize, ncol(m)))
    means = unlist(blockApply(m, function(x) mean(x, na.rm=na.rm), grid=grid, BPPARAM=BPPARAM))
    if(na.rm) {
        ns = unlist(blockApply(m, function(x) sum(!is.na(x)), grid=grid, BPPARAM=BPPARAM))
    } else {
        ns = length(m)
    }
    # weighted mean
    return(sum(means * (ns/sum(ns))))
}

block_rowMeans <- function(m, blockSize=10000, BPPARAM=SerialParam()) {
    grid = RegularArrayGrid(refdim=dim(m), spacings=c(blockSize, ncol(m)))
    return(unlist(blockApply(m, rowMeans, grid=grid, BPPARAM=BPPARAM)))
}

block_rowSums <- function(m, blockSize=10000, BPPARAM=SerialParam()) {
    grid = RegularArrayGrid(refdim=dim(m), spacings=c(blockSize, ncol(m)))
    return(unlist(blockApply(m, rowSums, grid=grid, BPPARAM=BPPARAM)))
}

block_rowVars <- function(m, blockSize=10000, BPPARAM=SerialParam()) {
    grid = RegularArrayGrid(refdim=dim(m), spacings=c(blockSize, ncol(m)))
    return(unlist(blockApply(m, rowVars, grid=grid, BPPARAM=BPPARAM)))
}

### some utility functions
## variance of each row
rowVars <- function(x) {
    rowMeans = DelayedArray::rowMeans
    n0 <- ncol(x)
    EX <- rowMeans(x, na.rm=TRUE)
    EX2 <- rowMeans(x^2, na.rm=TRUE)
    vv <- (EX2-EX^2) * n0 / (n0-1)
    vv
}

### generate negative binomial rv, given mu and phi (over-dispersion)
rnegbinom <- function (n, mu =1, phi=0.01){
    rpois(n, rgamma(n, shape=1/phi,scale=mu*phi))
}


###############################################
## Smoothing function. Smooth by chr
###############################################
smooth.chr <- function(x, ws, allchr, allpos, method=c("avg", "sum")) {
    method <- match.arg(method)
    if(method == "avg")
        flag = 1
    else
      flag = 0

    ## remove NA's
    ix = !is.na(x)
    x2 = x[ix]
    allchr2 = allchr[ix]
    allpos2 = allpos[ix]

    n0=length(x2)
    idx=split(1:n0, allchr2)
    res2=rep(0, n0)
    for(i in seq(along=idx)) {
        res2[idx[[i]]]=.Call("windowFilter", x2[idx[[i]]], as.integer(allpos2[idx[[i]]]), as.integer(ws), as.integer(flag))
    }

    if(sum(ix) > 0) {
        res = rep(NA, length(x))
        res[ix] = res2
    } else res = res2

    return(res)
}

