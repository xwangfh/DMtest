#' Compute p values based on DMC/DVC and joint tests.
#'
#' This function loads beta matrix, covariate matrix and design formula and
#' it computes p values based on DMC/DVC and joint tests.
#'
#' @param beta beta value matrix
#' @param covariate covariate matrix
#' @param npermut number of permutation
#' @param permut.seed random seed used for permutation
#' @param corenumber number of cores used for joint-tests
#' @return A data frame out for p-values.
#' DMCP: p-value from DMC test.
#' DVCP p-value from DVC test.
#' Joint1P: joint test for DMVC (test for
#' hypermethylation and increased variance in cancer).
#' Joint2P: joint test for DMVC (test for differential methylation
#' in both direction and increased variance in cancer).
#' @export
dmvc=function(beta=beta,covariate=covariate,npermut=100,permut.seed=100,corenumber=1)
{
  #preprocess data
  # check if covariate matrix has group variable
  if (sum(colnames(covariate)=="group")==0)
  {
    stop("covriate matrix should include group variable")
  }
  
  #make group to be the first covariate (or else it could be a problem for rare cases where some factors have one level and get removed for some CpGs due to a lot of missing data)
  tmp=colnames(covariate)
  tmp=tmp[tmp!="group"]
  tmp=c("group",tmp)
  covariate=covariate[,tmp]
  
  if (any(!covariate$group %in% c(0,1,NA)))
  {
    stop("group variable should be coded as 0 and 1")
  }
  #check if samples in beta match samples in covariate
  if ((!all(colnames(beta) %in% rownames(covariate))) | (! all(rownames(covariate) %in% colnames(beta))))
  {
    warning("There are extra samples in beta or covariate matrices")
  }
  allsamples=intersect(colnames(beta),rownames(covariate))
  if (length(allsamples)==0)
  {
    stop("There are no common samples in beta and covariate matrices")
  }
  idx=match(allsamples,colnames(beta))
  beta=beta[,idx]
  idx=match(allsamples,rownames(covariate))
  covariate=covariate[idx,]

  #remove NAs in covariate
  idx1=stats::complete.cases(covariate)
  if (sum(!idx1)>0)
  {
    warning(paste0("There are missing data in covariate matrix, ",sum(!idx1)," samples were removed"))
    beta=beta[,idx1]
    covariate=covariate[idx1,]
  }
  out=data.frame(Mean_normal=rep(NA,nrow(beta)),Mean_tumor=NA,Mean_all=NA,SD_normal=NA,SD_tumor=NA,
                 SD_all=NA,DMCP=rep(NA,nrow(beta)),DVCP=NA,Joint1P=NA,Joint2P=NA,LRT1=NA,LRT2=NA,pho=NA)
  rownames(out)=rownames(beta)
  beta[beta <= 0]=0.01
  beta[beta >= 1]=0.99
  out$Mean_normal=rowMeans(beta[,covariate$group==0],na.rm=T)
  out$Mean_tumor=rowMeans(beta[,covariate$group==1],na.rm=T)
  out$Mean_all=rowMeans(beta,na.rm=T)
  out$SD_normal=matrixStats::rowSds(as.matrix(beta[,covariate$group==0]),na.rm=T)
  out$SD_tumor=matrixStats::rowSds(as.matrix(beta[,covariate$group==1]),na.rm=T)
  out$SD_all=matrixStats::rowSds(as.matrix(beta),na.rm=T)
  #M values
  dat=log(beta/(1-beta),base=2)
  writeLines (paste0("Analyze ",nrow(beta)," CpGs across ",ncol(beta)," samples at ",as.character(Sys.time())))
  writeLines(paste0("Start to compute DMCP pvalues at ",as.character(Sys.time()),"..."))
  fit2all=fit.model.probe(mat = as.matrix(dat),colData = covariate)
  out$DMCP=fit2all$statistics$pval
  writeLines(paste0("Start to compute DVCP pvalues at ",as.character(Sys.time()),"..."))
  dat1=as.matrix(dat)
  dat1[,covariate$group==1] <- abs(dat1[,covariate$group==1]-matrixStats::rowMedians(dat1[,covariate$group==1],na.rm=T))
  dat1[,covariate$group==0] <- abs(dat1[,covariate$group==0]-matrixStats::rowMedians(dat1[,covariate$group==0],na.rm=T))
  fit1all=fit.model.probe(mat = dat1,colData = covariate)
  out$DVCP=fit1all$statistics$pval
  
  design=as.formula("~.")
  designmatrix=model.matrix(design,data=covariate)
  idxgroup=which(colnames(designmatrix)=="group")
  writeLines(paste0("Start to compute joint-test pvalues at ",as.character(Sys.time()),"..."))
 
  constrainMax <- function(coef.gg,cov.gg){
    dd <-  1/(cov.gg[1,1]*cov.gg[2,2]-cov.gg[1,2]^2)
    v11 <- cov.gg[2,2]
    v12 <- -cov.gg[1,2]
    v21 <- -cov.gg[2,1]
    v22 <- cov.gg[1,1]
    
    a1 <- v11
    b1 <- -2*v11*coef.gg[1]-coef.gg[2]*(v12+v21)
    c1 <- v11*coef.gg[1]^2+v22*coef.gg[2]^2+(v12+v21)*coef.gg[1]*coef.gg[2]
    
    a2 <- v22
    b2 <- -2*v22*coef.gg[2]-coef.gg[1]*(v12+v21)
    c2 <- v22*coef.gg[2]^2+v11*coef.gg[1]^2+(v12+v21)*coef.gg[1]*coef.gg[2]
    
    mlrt1 <- dd*ifelse(b1<0, c1-a1*((-b1)/(2*a1))^2,c1)
    mlrt2 <- dd*ifelse(b2<0, c2-a2*((-b2)/(2*a2))^2,c2)
    
    mlrt3 <- dd*(c2-a2*((-b2)/(2*a2))^2)
    
    return(list(lrt1=min(mlrt1,mlrt2),lrt2=mlrt3))
  }
  
  #use to permute y
  set.seed(permut.seed)
  #permutate index used for probes without missing data
  permutmat=matrix(NA,nrow=ncol(dat),ncol=npermut)
  for (j in 1:npermut)
  {
    # permutmat[,j]=sample(1:ncol(dat))
    permutmat[designmatrix[,idxgroup]==1,j]=sample((1:ncol(dat))[designmatrix[,idxgroup]==1],replace=T)
    permutmat[designmatrix[,idxgroup]==0,j]=sample((1:ncol(dat))[designmatrix[,idxgroup]==0],replace=T) 
  }
  #used for probes with missing data
  permutmatindx=function(y,npermut,permut.seed,designmatrixmis)
  {
    set.seed(permut.seed)
    permutmat=matrix(NA,nrow=length(y),ncol=npermut)
    for (j in 1:npermut)
    {
      #permutmat[,j]=sample(1:length(y))
      permutmat[designmatrixmis[,idxgroup]==1,j]=sample((1:length(y))[designmatrixmis[,idxgroup]==1],replace=T)
      permutmat[designmatrixmis[,idxgroup]==0,j]=sample((1:length(y))[designmatrixmis[,idxgroup]==0],replace=T)
    }
    return(permutmat)
  }
  #the permutated matrix
  yy1 <- matrix(NA,nrow=ncol(dat),ncol=npermut)
  #outbeta <- matrix(0,npermut,2)
  #iblock is the index of CpGs to work on at the same time
  compute_pho=function(dat,covariate,iblock=c(1,2),permutmat)
  {
    y0=dat[iblock,]
    y00=sapply(1:length(iblock), function(x) y0[x,][permutmat])
    npermut=ncol(permutmat)
    y10=matrix(y00,ncol=length(iblock)*npermut)
    yy10=matrix(NA,nrow=nrow(y10),ncol=ncol(y10))
    yy10[covariate$group==1,]=abs(sweep(y10[covariate$group==1,],2,matrixStats::colMedians(y10[covariate$group==1,],na.rm=T)))
    yy10[covariate$group==0,]=abs(sweep(y10[covariate$group==0,],2,matrixStats::colMedians(y10[covariate$group==0,],na.rm=T)))
    outbeta0 <- matrix(0,npermut*length(iblock),2)
    outbeta0[,1] <- matrixStats::colMeans2(yy10[covariate$group==1,],na.rm=T)- matrixStats::colMeans2(yy10[covariate$group==0,],na.rm=T)
    outbeta0[,2] <- matrixStats::colMeans2(y10[covariate$group==1,],na.rm=T)- matrixStats::colMeans2(y10[covariate$group==0,],na.rm=T)
    pho0=sapply(1:length(iblock),function(x) {
      segidx=seq((1+(x-1)*npermut),x*npermut)
      cor(outbeta0[segidx,1],outbeta0[segidx,2])})
    return(pho0)
  }
  
  pho_all=rep(NA,nrow(dat)) #store all pho
  #max nblock for 4G memory
  #4000*1e6/ncol(dat)/npermut/8
  nblock=100 #number of CpGs together for permutation
  nprobes=nrow(dat)
  if (nprobes<nblock) nblock=nprobes
  n.chunk=floor(nprobes/nblock)
  ncount <- rep(nblock,n.chunk)
  if (nprobes>sum(ncount))  ncount[n.chunk] <- ncount[n.chunk]+ nprobes-sum(ncount)
  chunks <- c(0,cumsum(ncount))
  
  if (corenumber>1) #parallel version for permutation
  {
    numberofcores=floor(parallel::detectCores()-1) #don't use all the cores
    if (corenumber>numberofcores) corenumber=numberofcores
    if (corenumber>1)
    {
      cl <- parallel::makeCluster(corenumber)
      doParallel::registerDoParallel(cl)
      if (corenumber>1) writeLines(paste0("Parallel version for permutation, total number of cores available:",numberofcores,", used cores:",corenumber))
      # library(bigstatsr)
      # pho_all=bigstatsr::FBM(nrow(dat),1)
      pho_all=foreach::foreach(r=1:n.chunk, .combine='c') %dopar% {
        iblocks=seq(chunks[r]+1,chunks[r+1])
        compute_pho(dat=dat,covariate=covariate,iblock=iblocks,permutmat=permutmat)
      }
      # writeLines(paste0("Permutation finishes ",as.character(Sys.time())))
      parallel::stopCluster(cl)
    }
  }
  if (corenumber<=1) #sequential version for permutation
  {
    for (r in 1:n.chunk)
    {
      iblocks=seq(chunks[r]+1,chunks[r+1])
      pho_all[iblocks]=compute_pho(dat=dat,covariate=covariate,iblock=iblocks,permutmat=permutmat)
    }
    # writeLines(paste0("Permutation part finishes ",as.character(Sys.time())))
  }
  tmp=sapply(1:nrow(dat),function(i){
    # for (i in 1:nrow(beta))
    # {
      out1=rep(NA,5)
      ## 2 df test estimating covariance ##
      pho=pho_all[i]
      cov.gg <- matrix(0,2,2)
      cov.gg[1,2] <-  fit1all$statistics$sd[i]*fit2all$statistics$sd[i]*pho #summary(fit1)$coef[2,2]*summary(fit2)$coef[2,2]*pho
      cov.gg[2,1] <-  cov.gg[1,2]
      cov.gg[1,1] <-  (fit1all$statistics$sd[i])^2
      cov.gg[2,2] <-  (fit2all$statistics$sd[i])^2
      coef.gg <- c(fit1all$statistics$coef[i],fit2all$statistics$coef[i])

      LRTtmp1=drop(crossprod(coef.gg,solve(cov.gg)) %*% coef.gg)
      constrainres=constrainMax(coef.gg,cov.gg)

      if (coef.gg[1]>=0 & coef.gg[2]>=0) LRT1 <- LRTtmp1 else LRT1<- LRTtmp1 - constrainres$lrt1
      #pho <- cov.gg[1,2]/sqrt(cov.gg[1,1]*cov.gg[2,2])
      if (is.numeric(LRT1))
      {
        out1[1] <- 0.5*pchisq(LRT1,df=1,lower.tail=F)+0.5*(1-(1/pi)*acos(pho))*pchisq(LRT1,df=2,lower.tail=F)
        out1[3]=LRT1
      }

      if (coef.gg[1]>=0) LRT2 <- LRTtmp1 else LRT2<- LRTtmp1 - constrainres$lrt2

      if (is.numeric(LRT2))
      {
        out1[2] <- 0.5*pchisq(LRT2,df=1,lower.tail=F)+0.5*pchisq(LRT2,df=2,lower.tail=F)
        out1[4]=LRT2
      }
      out1[5]=pho
      return(out1)
    }
    )
  
  tmp=t(tmp)
  out$Joint1P[1:nrow(tmp)]=tmp[,1]
  out$Joint2P[1:nrow(tmp)]=tmp[,2]
  out$LRT1[1:nrow(tmp)]=tmp[,3]
  out$LRT2[1:nrow(tmp)]=tmp[,4]
  out$pho[1:nrow(tmp)]=tmp[,5]
  writeLines(paste0("The result is ready at ",as.character(Sys.time())))
  return(out)
}

#test all CpGs
# Cancer.type <- "TCGA-COAD"
# RData.file <- paste("/fh/fast/dai_j/CancerGenomics/prostate_methylation/data/TCGA/",Cancer.type, "_", "methy.RData", sep='') # RData file for saving the TCGA data
# load(RData.file)
# 
# sampleids=c(colnames(normal_methy),colnames(tumor_methy_all))
# group <- c(rep(0,ncol(normal_methy)),rep(1,ncol(tumor_methy_all)))
# idx=match(sampleids,rownames(clinical))
# covariate=clinical[idx,]
# covariate=cbind(group,covariate)
# for (i in 1:ncol(covariate))
# {
#   if (is.character(covariate[,i])) covariate[,i]=as.factor(covariate[,i])
# }
# #designmatrix=model.matrix(~.,data=covariate)
# beta=cbind(normal_methy,tumor_methy_all)
# colnames(beta)=rownames(covariate)
