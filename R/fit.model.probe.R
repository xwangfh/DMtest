#' Fit linear model
#'
#' This function loads beta matrix, covariate matrix and
#' it computes statistics and fitted values.
#'
#' @param mat m value matrix
#' @param colData covariate matrix
#' @param design formula
#' @return A list including a matrix stats for test statistics and
#' a matrix fitted for fitted values
#' @export
fit.model.probe <- function(mat,colData){
  # Fit model for each probe from data
  # Call function "calculate.stats" to calculate coef & sd & residual values
  # Fast version using QR factorization
  # Employ the way of how glm handles NA using na.omit() by default
  # If our data have NAs, the probes with NAs will be calculated seperately with NA-omitted-data
  design=as.formula("~.")
  calculate.stats <- function(mat, colData){
    # Calculate coef & sd & residual values for selected probes from data
    # Only used by function "fit.model.probe"

    design.sub <- colData[rownames(colData) %in% colnames(mat),]
    #remove factors with 1 level
    idx_1lev=1:ncol(design.sub)
    idx_1lev=idx_1lev[sapply(idx_1lev,function(j) {length(unique(design.sub[,j]))==1})]

    #factors to be removed (0 means nothing to remove)
    removed.factors=0
    if (length(idx_1lev)>0)
    {
      design.sub=design.sub[,-idx_1lev]
      removed.factors=paste0(idx_1lev,collapse = ";")
    }
    design.sub <- model.matrix(design, data = design.sub)

    # Calculate coef
    M <- design.sub[,-2, drop=FALSE]
    M.qr <- qr(M)
    M.qr.Q <- qr.Q(M.qr)
    S <- diag(nrow(M)) - tcrossprod(M.qr.Q)
    V <- matrix(design.sub[,2], ncol=1)
    SV <- S %*% V
    coef <- (mat %*% crossprod(diag(nrow(M)) - tcrossprod(M.qr.Q), design.sub[,2])) / diag(crossprod(V,SV))

    # Calculate residuals
    qr.X <- qr(design.sub)
    qr.X.Q <- qr.Q(qr.X)
    resids <- tcrossprod(diag(nrow(design.sub)) - tcrossprod(qr.X.Q), mat)
    rownames(resids) <- colnames(mat)
    colnames(resids) <- rownames(mat)

    #Calculate fitted
    fitted =tcrossprod(tcrossprod(qr.X.Q),mat)
    rownames(fitted) <- colnames(mat)
    # Calculate SE
    qr.X.R <- qr.R(qr.X)
    SE <- sqrt(chol2inv(qr.X.R)[2,2] * (colSums(resids^2) / (ncol(mat) - M.qr$rank - 1)))

    # Calculate p-value
    pval= t(2*pt(abs(coef/SE), ncol(mat) - M.qr$rank - 1, lower.tail=FALSE))[1,]
    stats.sub <- as.matrix(cbind(coef, SE, coef / SE,pval,removed.factors))
    rownames(stats.sub) <- rownames(mat)
    colnames(stats.sub) <- c("coef", "sd", "zscore","pval","removedfactors")

    #
    # # Fill the fitted matrix/resids to full size of all samples. fitted of samples with NAs are NA.
    removed.samples <- rownames(colData)[!rownames(colData) %in% colnames(mat)]
    removed.mat <- matrix(NA, nrow = length(removed.samples), ncol = ncol(fitted))
    rownames(removed.mat) <- removed.samples
    fitted <- rbind(fitted, removed.mat)
    fitted <- fitted[rownames(colData), , drop = FALSE]
    
    resids <- rbind(resids, removed.mat)
    resids <- resids[rownames(colData), , drop = FALSE]

    #return(list(stats = stats.sub, residuals = resids))
    return(list(stats = stats.sub,fitted=fitted,design.sub=design.sub))
  }

  colData.factors=NULL
  for (i in 1:ncol(colData))
  {
    if (is.factor(colData[,i]))
    {
      colData.factors=c(colData.factors,i)
    }
  }

  NA.comb <- apply(mat, 1, function(x){paste(unname(which(is.na(x))), collapse = ';')})
  NA.comb.unique <- unique(NA.comb)

  stats.res <- lapply(NA.comb.unique, function(x){
    mat <- mat[NA.comb == x,, drop = FALSE]
    ##if(x != "") mat <- mat[, -as.numeric(strsplit(x, ';')[[1]]), drop = FALSE]
    ##change to, a probe have missing in multiple samples
    if(x != "") mat <- mat[, -as.numeric(unlist(strsplit(x, ';'))), drop = FALSE]
    stats.res.tmp <- calculate.stats(mat,colData)
    return(stats.res.tmp)
  })

  stats <- NULL

  stats$statistics <- do.call("rbind", lapply(stats.res, function(x){x$stats}))
  stats$statistics <- stats$statistics[rownames(mat),]
  stats$statistics=as.data.frame(stats$statistics)
  stats$statistics$coef=as.numeric(stats$statistics$coef)
  stats$statistics$sd=as.numeric(stats$statistics$sd)
  stats$statistics$zscore=as.numeric(stats$statistics$zscore)
  stats$statistics$pval=as.numeric(stats$statistics$pval)
  stats$fitted <- do.call("cbind", lapply(stats.res, function(x){x$fitted}))
  stats$fitted <- stats$fitted[,rownames(mat)]
  stats$fitted=as.data.frame(stats$fitted)
  stats$design=lapply(stats.res, function(x){x$design.sub})
  #designidx is used to access design matrix. stats$design[[designidx$designidx[i]]] is the design matrix for probe i (mat[i,])
  designidx=data.frame(matidx=1:nrow(mat),designidx=NA)
  designidx$designidx=match(NA.comb,NA.comb.unique)
  stats$designidx=designidx
  # stats$residuals <- do.call("cbind", lapply(stats.res, function(x){x$residuals}))
  # stats$residuals <- stats$residuals[,rownames(mat)]
  return(stats)
}

