readHeader<-function(dataset){
  data<-read.csv(dataset, header=F)
  header<-character()
  for (i in 1:ncol(data)){
    header[i]=toString(data[1,i])
  }
  return(header)
}
Read<-function(dataset,sep=","){
  data<-read.csv(dataset, header=TRUE, sep=sep)
  return(data)
}
Standardise<-function(data){
  ncol<-ncol(data)
  nrow<-nrow(data)
  #stdData<-matrix(nrow=nrow,ncol=ncol-1)
  stdData<-matrix(nrow=nrow,ncol=ncol)
  rownames(stdData)<-data[,1]
  #for (i in 2:ncol){
  for (i in 1:ncol){
    #stdData[,i-1]<-scale(data[i],center=TRUE, scale=TRUE)
    stdData[,i]<-scale(data[i],center=TRUE, scale=TRUE)
    
  }
  return(stdData)
}

geneSymbol<-function(mRNA){
  library(AnnotationFuncs,quietly=TRUE)
  library(org.Hs.eg.db,quietly=TRUE)
  #convert genesymbol to entrezID
  symbol<-character()
  for (i in 1:length(mRNA)){
    symbol[i]<-toString(translate(mRNA[i],org.Hs.egALIAS2EG)) #translate gene symboe to EntrezID.0 (Homo Sapien!)
  }
  return(symbol)
}
match<-function(miIndex,gene,targetall){ #find if gene is the target of the (miIndex)th  miRNA 
  #if (is.na(targetall[miIndex]))
  #  result=FALSE
  #else
{
  targetall[[miIndex]]<-targetall[[miIndex]][!is.na(targetall[[miIndex]])]
  for (i in 1:length(targetall[[miIndex]])){
    
    if (gene==(targetall[[miIndex]][i])){
      result=TRUE
      break
    }
    else
    {result=FALSE}
  }
}
return(result)
}
getEdge<-function(targetall,geneSymbol,num_miRNA,nrow,ncol){   
  #targetall: a list of miRNA targets
  #geneSymbol: a character array contains all gene symbols
  edge<-matrix(FALSE,nrow,ncol)
  for (i in 1:length(miRNA)){
    print(i)
    if (is.na(targetall)[i]==TRUE){
      next
    }
    else{
      for (j in 1:length(geneSymbol)){
        if (match(i,geneSymbol[j],targetall))
        {edge[i,j+num_miRNA]=TRUE;
         #          edge[j+35,i]=TRUE
        }
      }
    } 
  }
  return(edge)
}
queryTarget<-function(miRNA){
  library("org.Bt.eg.db",quietly=T)
  library("targetscan.Hs.eg.db",quietly=T)
  name=names(data)
  #query database for microRNA target
  targetall<-character()
  familyname<-character()
  for (i in 1:length(miRNA)){
    familyname<-try(toString(mget(miRNA[i], targetscan.Hs.egMIRBASE2FAMILY))) #get the miR family name
    if(class(familyname) == "try-error") {targetall[i]<-NA;next}
    targetall[i]<-mget(familyname, revmap(targetscan.Hs.egTARGETS)) #get the target for miR family(EntrezID)
  }
  
  return(targetall)
}
queryTarget2<-function(miRNA,database="targetscan"){
  library(RmiR.Hs.miRNA,quietly=T)
  target<-list()
  for (i in 1:length(miRNA)){
    temp<-dbGetQuery(RmiR.Hs.miRNA_dbconn(), paste("SELECT * FROM ",database," WHERE mature_miRNA='",miRNA[i],"'",sep="")) #string concatenation
    if (length(temp$gene_id)>1){
      target[[i]]<-as.character(temp[,2])}   #target list for a particular miRNA
    else
    {target[i]<-NA;next}
  }
  return(target)
}


queryTargetFile<-function(miRNA,mRNA,file){
  #the column name of file should be tagged as "mir" and "gene"
  data<-Read(file)
  mir=as.character(data$mir)
#   mir<-paste("hsa-",sep="",mir);mir<-sub('r','R',mir)
  gene=as.character(data$gene)
#   symbol<-geneSymbol(gene)
  sum=0
  rep<-replicate(length(miRNA),mir)
  edge=matrix(F,length(miRNA)+length(mRNA),length(miRNA)+length(mRNA))
  for (i in 1:length(mir)){
#     print(i)
    if (length(which(rep[i,]==miRNA)>0)){
      match1<-which(rep[i,]==miRNA,arr.ind=T)
      #gene: gene[i] #mirna: miRNA[match]
      rep2<-replicate(length(mRNA),gene[i])
      match2<-which(rep2==mRNA,arr.ind=T)
      edge[match1,match2+length(miRNA)]=T
    }
  }
return(edge)
}


queryPPIFile<-function(miRNA,mRNA,file){
  #the column name of file should be tagged as "mir" and "gene"
  data<-read.csv(file,as.is=T)
  gene1=data$gene1
  gene2=data$gene2
  sum=0
  rep<-replicate(length(mRNA),gene1)
  edge=matrix(F,length(miRNA)+length(mRNA),length(miRNA)+length(mRNA))
  for (i in 1:length(mRNA)){
    if (length(which(rep[i,]==mRNA)>0)){
      match1<-which(rep[i,]==mRNA,arr.ind=T)
      #gene: gene[i] #mirna: miRNA[match]
      rep2<-replicate(length(mRNA),gene2[i])
      match2<-which(rep2==mRNA,arr.ind=T)
      edge[length(miRNA)+match1,match2+length(miRNA)]=T
    }
  }
  return(edge)
}



CausaleffectsPC<-function(miRNA,mRNA,stdData,alp, nmiRs,fixedEdge=NULL,solve.confl=FALSE,skel.method="original"){
  
  #Learning PCDAG
  suffStat<-list(C=cor(stdData), n=nrow(stdData))  
  
  print('Begin PC')
  pc.fit<-pc(suffStat, indepTest=gaussCItest, p=ncol(stdData), alpha=alp,fixedEdge=fixedEdge,skel.method="original",solve.confl=solve.confl)               
  print("PC finished")
  
  nmRs<-ncol(stdData)-nmiRs
  result<-matrix(nrow=nmRs, ncol=nmiRs)
  #for (l in 1:nmiRs){
  l=1
  while(l<=nmiRs){
    #Inferring causal effects
    caef<-try(idaFast(l,(nmiRs+1):ncol(stdData),cov(stdData), pc.fit@graph ))
    if(class(caef) == "try-error") {result=NA;print("singular error");break}
    else{
      #min of absolute values.
      caef1<-matrix(nrow=nmRs,ncol=1)
      for (k in 1:nmRs){
        caefabs<-abs(caef)
        index<-which(caefabs==min(caefabs[k,]), arr.ind=TRUE)
        pos<-index[1,2]
        caef1[k,]<-caef[k,pos]
      }
      result[,l]<-caef1
      l=l+1
    }
  }
  rownames(result)<-mRNA
  colnames(result)<-miRNA
  print('Causaul Finish!')
  return(result)
}
CausaleffectsWeijia<-function(miRNA,mRNA,stdData,alp, nmiRs,fixedEdge=NULL,solve.confl=FALSE,skel.method=skel.method){

  #Learning PCDAG
  suffStat<-list(C=cor(stdData), n=nrow(stdData))  
  
  print('Begin PC')
  pc.fit<-modifiedpcWeijia(suffStat, indepTest=gaussCItest, p=ncol(stdData), alpha=alp,fixedEdge=fixedEdge,skel.method=skel.method,solve.confl=solve.confl)               
  print("PC finished")
  
  nmRs<-ncol(stdData)-nmiRs
  result<-matrix(nrow=nmRs, ncol=nmiRs)
  #for (l in 1:nmiRs){
  l=1
  while(l<=nmiRs){
    #Inferring causal effects
    caef<-try(idaFast(l,(nmiRs+1):ncol(stdData),cov(stdData), pc.fit@graph ))
    print(l)
    if(class(caef) == "try-error") {result=NA;print("singular error");break}
    else{
      #min of absolute values.
      caef1<-matrix(nrow=nmRs,ncol=1)
      for (k in 1:nmRs){
        caefabs<-abs(caef)
        index<-which(caefabs==min(caefabs[k,]), arr.ind=TRUE)
        pos<-index[1,2]
        caef1[k,]<-caef[k,pos]
      }
      result[,l]<-caef1
      l=l+1
    }
  }
  rownames(result)<-mRNA
  colnames(result)<-miRNA
  print('Causaul Finish!')
  return(result)
}
CausaleffectsThuc<-function(miRNA,mRNA,stdData,alp, nmiRs,fixedEdge=NULL,u2pd = c("relaxed", "rand", "retry"),solve.confl = FALSE,skel.method="original"){
  
  #Learning PCDAG
  suffStat<-list(C=cor(stdData), n=nrow(stdData))  
  
  print('Begin PC')
  pc.fit<-modifiedpcThuc(suffStat, indepTest=gaussCItest, p=ncol(stdData), alpha=alp,fixedEdge=fixedEdge,skel.method=skel.method,u2pd = u2pd,solve.confl=solve.confl)               
  print("PC finished")
  
  nmRs<-ncol(stdData)-nmiRs
  result<-matrix(nrow=nmRs, ncol=nmiRs)
  #for (l in 1:nmiRs){
  l=1
  while(l<=nmiRs){
    #Inferring causal effects
    caef<-try(idaFast2(l,(nmiRs+1):ncol(stdData),cov(stdData), pc.fit@graph ))
    if(class(caef) == "try-error") {result=NA;print("singular error");break}
    else{
      #min of absolute values.
      caef1<-matrix(nrow=nmRs,ncol=1)
      for (k in 1:nmRs){
        caefabs<-abs(caef)
        index<-which(caefabs==min(caefabs[k,]), arr.ind=TRUE)
        pos<-index[1,2]
        caef1[k,]<-caef[k,pos]
      }
      result[,l]<-caef1
      l=l+1
    }
  }
  rownames(result)<-mRNA
  colnames(result)<-miRNA
  print('Causaul Finish!')
  return(result)
}


skeletonWeijia <- function(suffStat, indepTest, alpha, labels, p,
                     method = c("stable", "original", "stable.fast"), m.max = Inf,
                     fixedGaps = NULL, fixedEdges = NULL,
                     NAdelete = TRUE, verbose = FALSE){
  cl <- match.call()
  if(!missing(p)) stopifnot(is.numeric(p), length(p <- as.integer(p)) == 1, p >= 2)
  if(missing(labels)) {
    if(missing(p)) stop("need to specify 'labels' or 'p'")
    labels <- as.character(seq_len(p))
  } else { ## use labels ==> p  from it
    stopifnot(is.character(labels))
    if(missing(p)) {
      p <- length(labels)
    } else if(p != length(labels))
      stop("'p' is not needed when 'labels' is specified, and must match length(labels)")
    else
      message("No need to specify 'p', when 'labels' is given")
  }
  seq_p <- seq_len(p)
  method <- match.arg(method)
  
  ## G := !fixedGaps, i.e. G[i,j] is true  iff  i--j  will be investigated
  if (is.null(fixedGaps)) {
    G <- matrix(TRUE, nrow = p, ncol = p)
    diag(G) <- FALSE
  }
  else if (!identical(dim(fixedGaps), c(p, p)))
    stop("Dimensions of the dataset and fixedGaps do not agree.")
  else if (!identical(fixedGaps, t(fixedGaps)) )
    stop("fixedGaps must be symmetric")
  else
    G <- !fixedGaps
  
  if (any(is.null(fixedEdges))) { ## MM: could be sparse
    fixedEdges <- matrix(rep(FALSE, p * p), nrow = p, ncol = p)
  }
  else if (!identical(dim(fixedEdges), c(p, p)))
    stop("Dimensions of the dataset and fixedEdges do not agree.")

  
  if (method == "stable.fast") {
    ## Do calculation in C++...
    if (identical(indepTest, gaussCItest))
      indepTestName <- "gauss"
    else
      indepTestName <- "rfun"
    options <- list(
      verbose = as.integer(verbose), 
      m.max = as.integer(ifelse(is.infinite(m.max), p, m.max)),
      NAdelete = NAdelete)
    res <- .Call("estimateSkeleton", G, suffStat, indepTestName, indepTest, alpha, fixedEdges, options);
    G <- res$amat
    # sepset <- res$sepset
    sepset <- lapply(seq_p, function(i) c(
      lapply(res$sepset[[i]], function(v) if(identical(v, as.integer(-1))) NULL else v),
      vector("list", p - length(res$sepset[[i]])))) # TODO change convention: make sepset triangular
    pMax <- res$pMax
    n.edgetests <- res$n.edgetests
    ord <- length(n.edgetests) - 1
  }
  else {
    ## Original R version
    G1=G-fixedEdges; G1=G1!=0
    pval <- NULL
    sepset <- lapply(seq_p, function(.) vector("list",p))# a list of lists [p x p]
    ## save maximal p value
    pMax <- matrix(-Inf, nrow = p, ncol = p)
    diag(pMax) <- 1
    done <- FALSE
    ord <- 0
    n.edgetests <- numeric(1)# final length = max { ord}
    while (!done && any(G1) && ord <= m.max) {
      n.edgetests[ord1 <- ord+1L] <- 0
      done <- TRUE
      ind <- which(G1, arr.ind = TRUE)
      ## For comparison with C++ sort according to first row
      ind <- ind[order(ind[, 1]), ]
      remainingEdgeTests <- nrow(ind)
      cat("Order=", ord, "; remaining edges:", remainingEdgeTests,"\n",sep="")
      if(method == "stable") {
        ## Order-independent version: Compute the adjacency sets for any vertex
        ## Then don't update when edges are deleted
        G1.l <- split(G1, gl(p,p))
      }
      for (i in 1:remainingEdgeTests) { 
        if (verbose && i%%100 == 0) cat("|i=", i, "|iMax=", nrow(ind), "\n")
        x <- ind[i, 1]
        y <- ind[i, 2]
        if (G1[y, x] && !fixedEdges[y, x]) {
          nbrsBool <- if(method == "stable") G1.l[[x]] else G1[,x]
          nbrsBool[y] <- FALSE
          nbrs <- seq_p[nbrsBool]
          length_nbrs <- length(nbrs)
          if (length_nbrs >= ord) {
            if (length_nbrs > ord)
              done <- FALSE
            S <- seq_len(ord)
            repeat { ## condition w.r.to all  nbrs[S] of size 'ord'
              n.edgetests[ord1] <- n.edgetests[ord1] + 1
              pval <- indepTest(x, y, nbrs[S], suffStat)
              if (verbose)
                cat("x=", x, " y=", y, " S=", nbrs[S], ": pval =", pval, "\n")
              if(is.na(pval))
                pval <- as.numeric(NAdelete) ## = if(NAdelete) 1 else 0
              if (pMax[x, y] < pval)
                pMax[x, y] <- pval
              if(pval >= alpha) { # independent
                G1[y, x]<-G1[x,y] <- FALSE
                sepset[[x]][[y]] <- nbrs[S]
                break
              }
              else {
                nextSet <- getNextSet(length_nbrs, ord, S)
                if (nextSet$wasLast)
                  break
                S <- nextSet$nextSet
              }
            } ## repeat
          }
        }
      }# for( i )
      ord <- ord + 1
    } ## while()
    for (i in 1:(p - 1)) {
      for (j in 2:p)
        pMax[i, j] <- pMax[j, i] <- max(pMax[i, j], pMax[j,i])
    }
  }
  
  ## transform matrix to graph object :
  Gobject <-
    if (sum(G1) == 0) {
      new("graphNEL", nodes = labels)
    } else {
      colnames(G1) <- rownames(G1) <- labels
      as(G1+fixedEdges,"graphNEL")
    }
#   return(G)
  # final object
  new("pcAlgo", graph = Gobject, call = cl, n = integer(0),
      max.ord = as.integer(ord - 1), n.edgetests = n.edgetests,
      sepset = sepset, pMax = pMax, zMin = matrix(NA, 1, 1))
  
}
skeletonDirectionsWeijia <- function(suffStat, indepTest, alpha, labels, p,
                           method = c("stable", "original", "stable.fast"), m.max = Inf,
                           fixedGaps = NULL, fixedEdges = NULL,
                           NAdelete = TRUE, verbose = FALSE,fixedDirections = NULL){
  cl <- match.call()
  if(!missing(p)) stopifnot(is.numeric(p), length(p <- as.integer(p)) == 1, p >= 2)
  if(missing(labels)) {
    if(missing(p)) stop("need to specify 'labels' or 'p'")
    labels <- as.character(seq_len(p))
  } else { ## use labels ==> p  from it
    stopifnot(is.character(labels))
    if(missing(p)) {
      p <- length(labels)
    } else if(p != length(labels))
      stop("'p' is not needed when 'labels' is specified, and must match length(labels)")
    else
      message("No need to specify 'p', when 'labels' is given")
  }
  seq_p <- seq_len(p)
  method <- match.arg(method)
  
  ## G := !fixedGaps, i.e. G[i,j] is true  iff  i--j  will be investigated
  if (is.null(fixedGaps)) {
    G <- matrix(TRUE, nrow = p, ncol = p)
    diag(G) <- FALSE
  }
  else if (!identical(dim(fixedGaps), c(p, p)))
    stop("Dimensions of the dataset and fixedGaps do not agree.")
  else if (!identical(fixedGaps, t(fixedGaps)) )
    stop("fixedGaps must be symmetric")
  else
    G <- !fixedGaps
  
  if (any(is.null(fixedEdges))) { ## MM: could be sparse
    fixedEdges <- matrix(rep(FALSE, p * p), nrow = p, ncol = p)
  }
  else if (!identical(dim(fixedEdges), c(p, p)))
    stop("Dimensions of the dataset and fixedEdges do not agree.")
  
  if (method == "stable.fast") {
    ## Do calculation in C++...
    if (identical(indepTest, gaussCItest))
      indepTestName <- "gauss"
    else
      indepTestName <- "rfun"
    options <- list(
      verbose = as.integer(verbose), 
      m.max = as.integer(ifelse(is.infinite(m.max), p, m.max)),
      NAdelete = NAdelete)
    res <- .Call("estimateSkeleton", G, suffStat, indepTestName, indepTest, alpha, fixedEdges, options);
    G <- res$amat
    # sepset <- res$sepset
    sepset <- lapply(seq_p, function(i) c(
      lapply(res$sepset[[i]], function(v) if(identical(v, as.integer(-1))) NULL else v),
      vector("list", p - length(res$sepset[[i]])))) # TODO change convention: make sepset triangular
    pMax <- res$pMax
    n.edgetests <- res$n.edgetests
    ord <- length(n.edgetests) - 1
  }
  else {
    ## Original R version
    G1=G-fixedEdges; G1=G1!=0
    pval <- NULL
    sepset <- lapply(seq_p, function(.) vector("list",p))# a list of lists [p x p]
    ## save maximal p value
    pMax <- matrix(-Inf, nrow = p, ncol = p)
    diag(pMax) <- 1
    done <- FALSE
    ord <- 0
    n.edgetests <- numeric(1)# final length = max { ord}
    while (!done && any(G1) && ord <= m.max) {
      n.edgetests[ord1 <- ord+1L] <- 0
      done <- TRUE
      ind <- which(G1, arr.ind = TRUE)
      ## For comparison with C++ sort according to first row
      ind <- ind[order(ind[, 1]), ]
      remainingEdgeTests <- nrow(ind)
      cat("Order=", ord, "; remaining edges:", remainingEdgeTests,"\n",sep="")
      if(method == "stable") {
        ## Order-independent version: Compute the adjacency sets for any vertex
        ## Then don't update when edges are deleted
        G1.l <- split(G1, gl(p,p))
      }
      for (i in 1:remainingEdgeTests) { 
        if (verbose && i%%100 == 0) cat("|i=", i, "|iMax=", nrow(ind), "\n")
        x <- ind[i, 1]
        y <- ind[i, 2]
        if (G1[y, x] && !fixedEdges[y, x]) {
          nbrsBool <- if(method == "stable") G1.l[[x]] else G1[,x]
          nbrsBool[y] <- FALSE
          nbrs <- seq_p[nbrsBool]
          length_nbrs <- length(nbrs)
          if (length_nbrs >= ord) {
            if (length_nbrs > ord)
              done <- FALSE
            S <- seq_len(ord)
            repeat { ## condition w.r.to all  nbrs[S] of size 'ord'
              n.edgetests[ord1] <- n.edgetests[ord1] + 1
              pval <- indepTest(x, y, nbrs[S], suffStat)
              if (verbose)
                cat("x=", x, " y=", y, " S=", nbrs[S], ": pval =", pval, "\n")
              if(is.na(pval))
                pval <- as.numeric(NAdelete) ## = if(NAdelete) 1 else 0
              if (pMax[x, y] < pval)
                pMax[x, y] <- pval
              if(pval >= alpha) { # independent
                G1[y, x]<-G1[x,y] <- FALSE
                sepset[[x]][[y]] <- nbrs[S]
                break
              }
              else {
                nextSet <- getNextSet(length_nbrs, ord, S)
                if (nextSet$wasLast)
                  break
                S <- nextSet$nextSet
              }
            } ## repeat
          }
        }
      }# for( i )
      ord <- ord + 1
    } ## while()
    for (i in 1:(p - 1)) {
      for (j in 2:p)
        pMax[i, j] <- pMax[j, i] <- max(pMax[i, j], pMax[j,i])
    }
  }
  ## transform matrix to graph object :
  Gobject <-
    if (sum(G1) == 0) {
      new("graphNEL", nodes = labels)
    } else {
      colnames(G1) <- rownames(G1) <- labels
      as(G1+fixedEdges-t(fixedDirections),"graphNEL")
    }
  #   return(G)
  # final object
  new("pcAlgo", graph = Gobject, call = cl, n = integer(0),
      max.ord = as.integer(ord - 1), n.edgetests = n.edgetests,
      sepset = sepset, pMax = pMax, zMin = matrix(NA, 1, 1))
  
}
idaFast2 <- function(x.pos, y.pos.set, mcov, graphEst)
{
  
  ## prepare adjMatrix and skeleton
  amat <- ad.g <- wgtMatrix.0(graphEst)   #wgtMatrix.0 modified 17/10/2014
  
  amat[which(amat!=0)] <- 1 ## i->j if amat[j,i]==1
  amatSkel <- amat + t(amat)
  amatSkel[amatSkel!=0] <- 1
  
  ## find unique and ambiguous parents of x
  wgt.est <- (ad.g != 0)
  tmp <- wgt.est-t(wgt.est)
  tmp[which(tmp<0)] <- 0
  wgt.unique <- tmp
  wgt.ambig <- wgt.est-wgt.unique
  pa1 <- which(wgt.unique[x.pos,]!=0)
  pa2 <- which(wgt.ambig[x.pos,]!=0)
  
  ## estimate beta
  if (length(pa2)==0) {
    print("lm.cov0")
    beta.tmp <- lm.cov(mcov,y.pos.set,c(x.pos,pa1)) ####
    beta.tmp[y.pos.set %in% pa1] <- 0
    beta.hat <- cbind(beta.tmp)
  } else {    ## at least one undirected parent
    ## no member of pa2
    pa2.f <- pa2
    pa2.t <- NA
    beta.hat <-
      if (!has.new.coll(amat,amatSkel,x.pos,pa1,pa2.t,pa2.f)) {
        print("lm.cov1")
        beta.tmp <- lm.cov(mcov,y.pos.set,c(x.pos,pa1)) ####
        beta.tmp[y.pos.set %in% pa1] <- 0
        cbind(beta.tmp)
      } # else NULL
    
    ## exactly one member of pa2
    for (i2 in seq_along(pa2)) {
      pa2.f <- pa2[-i2]
      pa2.t <- pa2[i2]
      if (!has.new.coll(amat,amatSkel,x.pos,pa1,pa2.t,pa2.f)) {
        print("lm.cov2")
        beta.tmp <- lm.cov(mcov,y.pos.set,c(x.pos,pa1,pa2.t)) ####
        beta.tmp[y.pos.set %in% c(pa1,pa2.t)] <- 0
        beta.hat <- cbind(beta.hat, beta.tmp)
      }
    } ## for (i2 in seq_along(pa2))
    
    ## higher order subsets of pa2
    if (length(pa2) > 1)
      for (i in 2:length(pa2)) {
        pa.tmp <- combn(pa2,i,simplify=TRUE)
        for (j in seq_len(ncol(pa.tmp))) {
          pa2.t <- pa.tmp[,j]
          pa2.f <- setdiff(pa2, pa2.t)
          if (!has.new.coll(amat,amatSkel,x.pos,pa1,pa2.t,pa2.f)) {
            beta.tmp <- lm.cov(mcov,y.pos.set,c(x.pos,pa1,pa2.t)) ####
            beta.tmp[y.pos.set %in% c(pa1,pa2.t)] <- 0
            beta.hat <- cbind(beta.hat, beta.tmp)
          }
        } ## for (j )
      } ## for (i )
    
  } ## if .. else length(pa2) > 0)
  
  ## MM: for now, maybe in the future get sensible column names:
  colnames(beta.hat) <- NULL
  if (nrow(beta.hat) > 0) rownames(beta.hat) <- as.character(y.pos.set)
  beta.hat
}
idaFast <- function(x.pos, y.pos.set, mcov, graphEst)
{
  ## Purpose: Estimate the causal effect of x on each element in the
  ## set y using the local method; graphEst and correlation matrix
  ## have to be precomputed; orient
  ## undirected edges at x in a way so that no new collider is
  ## introduced; if there is an undirected edge between x and y, both directions are considered;
  ## i.e., y might be partent of x in which case the effect is 0.
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - x.pos, y.pos: Column of x and y in d.mat
  ## - mcov: Covariance matrix that was used to estimate graphEst
  ## - graphEst: Fit of PC Algorithm (semidirected)
  ## ----------------------------------------------------------------------
  ## Value: list of causal values; one list element for each element of
  ## y.pos.set
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 7 Jan 2010, 11:18
  
  ## prepare adjMatrix and skeleton
  amat <- ad.g <- wgtMatrix.0(graphEst) #wgtMatrix.0 modified 17/10/2014
  
  amat[which(amat!=0)] <- 1 ## i->j if amat[j,i]==1
  amatSkel <- amat + t(amat)
  amatSkel[amatSkel!=0] <- 1
  
  ## find unique and ambiguous parents of x
  wgt.est <- (ad.g != 0)
  tmp <- wgt.est-t(wgt.est)
  tmp[which(tmp<0)] <- 0
  wgt.unique <- tmp
  wgt.ambig <- wgt.est-wgt.unique
  pa1 <- which(wgt.unique[x.pos,]!=0)
  pa2 <- which(wgt.ambig[x.pos,]!=0)
  
  ## estimate beta
  if (length(pa2)==0) {
    beta.tmp <- lm.cov(mcov,y.pos.set,c(x.pos,pa1)) ####
    beta.tmp[y.pos.set %in% pa1] <- 0
    beta.hat <- cbind(beta.tmp)
#     print("pa2==0")
  } else {    ## at least one undirected parent
    ## no member of pa2
#     print("at least one undirected parent no member of pa2")
    pa2.f <- pa2
    pa2.t <- NA
    beta.hat <-
      if (!has.new.coll(amat,amatSkel,x.pos,pa1,pa2.t,pa2.f)) {
#         print("if1")
        beta.tmp <- lm.cov(mcov,y.pos.set,c(x.pos,pa1)) ####
        beta.tmp[y.pos.set %in% pa1] <- 0
        cbind(beta.tmp)
      } # else NULL
    
    ## exactly one member of pa2
    for (i2 in seq_along(pa2)) {
#       print("## exactly one member of pa2")
      pa2.f <- pa2[-i2]
      pa2.t <- pa2[i2]
      if (!has.new.coll(amat,amatSkel,x.pos,pa1,pa2.t,pa2.f)) {
#         print("if2")
        beta.tmp <- lm.cov(mcov,y.pos.set,c(x.pos,pa1,pa2.t)) ####
        beta.tmp[y.pos.set %in% c(pa1,pa2.t)] <- 0
        beta.hat <- cbind(beta.hat, beta.tmp)
      }
    } ## for (i2 in seq_along(pa2))
    
    ## higher order subsets of pa2
    if (length(pa2) > 1)
#       print("    ## higher order subsets of pa2")
      for (i in 2:length(pa2)) {
        pa.tmp <- combn(pa2,i,simplify=TRUE)
        for (j in seq_len(ncol(pa.tmp))) {
          pa2.t <- pa.tmp[,j]
          pa2.f <- setdiff(pa2, pa2.t)
          if (!has.new.coll(amat,amatSkel,x.pos,pa1,pa2.t,pa2.f)) {
#             print("if3")
            beta.tmp <- lm.cov(mcov,y.pos.set,c(x.pos,pa1,pa2.t)) ####
            beta.tmp[y.pos.set %in% c(pa1,pa2.t)] <- 0
            beta.hat <- cbind(beta.hat, beta.tmp)
          }
        } ## for (j )
      } ## for (i )
    
  } ## if .. else length(pa2) > 0)
#   save(beta.hat,file="test.Rdata")
  ## MM: for now, maybe in the future get sensible column names:
  colnames(beta.hat) <- NULL
  if (nrow(beta.hat) > 0) rownames(beta.hat) <- as.character(y.pos.set)
  beta.hat
}


modifiedpcWeijia <- function(suffStat, indepTest, alpha, labels, p,
                       fixedGaps = NULL, fixedEdges = NULL, NAdelete = TRUE, m.max = Inf,
                       u2pd = c("relaxed", "rand", "retry"),
                       skel.method = c("stable", "original", "stable.fast"),
                       conservative = FALSE, maj.rule = FALSE,
                       solve.confl = FALSE, verbose = FALSE){  
  ## Initial Checks
  cl <- match.call()
  if(!missing(p)) stopifnot(is.numeric(p), length(p <- as.integer(p)) == 1, p >= 2)
  if(missing(labels)) {
    if(missing(p)) stop("need to specify 'labels' or 'p'")
    labels <- as.character(seq_len(p))
  } else { ## use labels ==> p  from it
    stopifnot(is.character(labels))
    if(missing(p)) {
      p <- length(labels)
    } else if(p != length(labels))
      stop("'p' is not needed when 'labels' is specified, and must match length(labels)")
    else
      message("No need to specify 'p', when 'labels' is given")
  }
  seq_p <- seq_len(p)
  
  u2pd <- match.arg(u2pd)
  skel.method <- match.arg(skel.method)
  if(u2pd != "relaxed") {
    if (conservative || maj.rule)
      stop("Conservative PC and majority rule PC can only be run with 'u2pd = relaxed'")
    
    if (solve.confl)
      stop("Versions of PC using lists for the orientation rules (and possibly bi-directed edges)\n can only be run with 'u2pd = relaxed'")
  }
  
  if (conservative && maj.rule) stop("Choose either conservative PC or majority rule PC!")
  
  ## Skeleton
  skel <- skeletonWeijia(suffStat, indepTest, alpha, labels=labels, method = skel.method,
                   fixedGaps=fixedGaps, fixedEdges=fixedEdges,
                   NAdelete=NAdelete, m.max=m.max, verbose=verbose)
  skel@call <- cl # so that makes it into result
  
  ## Orient edges
  if (!conservative && !maj.rule) {
    switch (u2pd,
            "rand" = udag2pdag(skel),
            "retry" = udag2pdagSpecial(skel)$pcObj,
            "relaxed" = udag2pdagRelaxed(skel, verbose=verbose, solve.confl=solve.confl))
  }
  else { ## u2pd "relaxed" : conservative _or_ maj.rule
    
    ## version.unf defined per default
    ## Tetrad CPC works with version.unf=c(2,1)
    ## see comment on pc.cons.intern for description of version.unf
    pc. <- pc.cons.intern(skel, suffStat, indepTest, alpha,
                          version.unf=c(2,1), maj.rule=maj.rule, verbose=verbose)
    udag2pdagRelaxed(pc.$sk, verbose=verbose,
                     unfVect=pc.$unfTripl, solve.confl=solve.confl)
  }
}
modifiedpcDirectionsWeijia <- function(suffStat, indepTest, alpha, labels, p,
                             fixedGaps = NULL, fixedEdges = NULL, fixedDirections = NULL, NAdelete = TRUE, m.max = Inf,
                             u2pd = c("relaxed", "rand", "retry"),
                             skel.method = c("stable", "original", "stable.fast"),
                             conservative = FALSE, maj.rule = FALSE,
                             solve.confl = FALSE, verbose = FALSE){  
  ## Initial Checks
  cl <- match.call()
  if(!missing(p)) stopifnot(is.numeric(p), length(p <- as.integer(p)) == 1, p >= 2)
  if(missing(labels)) {
    if(missing(p)) stop("need to specify 'labels' or 'p'")
    labels <- as.character(seq_len(p))
  } else { ## use labels ==> p  from it
    stopifnot(is.character(labels))
    if(missing(p)) {
      p <- length(labels)
    } else if(p != length(labels))
      stop("'p' is not needed when 'labels' is specified, and must match length(labels)")
    else
      message("No need to specify 'p', when 'labels' is given")
  }
  seq_p <- seq_len(p)
  
  u2pd <- match.arg(u2pd)
  skel.method <- match.arg(skel.method)
  if(u2pd != "relaxed") {
    if (conservative || maj.rule)
      stop("Conservative PC and majority rule PC can only be run with 'u2pd = relaxed'")
    
    if (solve.confl)
      stop("Versions of PC using lists for the orientation rules (and possibly bi-directed edges)\n can only be run with 'u2pd = relaxed'")
  }
  
  if (conservative && maj.rule) stop("Choose either conservative PC or majority rule PC!")
  
  ## Skeleton
  skel <- skeletonDirectionsWeijia(suffStat, indepTest, alpha, labels=labels, method = skel.method,
                         fixedGaps=fixedGaps, fixedEdges=fixedEdges, fixedDirections=fixedDirections,
                         NAdelete=NAdelete, m.max=m.max, verbose=verbose)
  skel@call <- cl # so that makes it into result
  
  ## Orient edges
  if (!conservative && !maj.rule) {
    switch (u2pd,
            "rand" = udag2pdag(skel),
            "retry" = udag2pdagSpecial(skel)$pcObj,
            "relaxed" = udag2pdagRelaxed(skel, verbose=verbose, solve.confl=solve.confl))
  }
  else { ## u2pd "relaxed" : conservative _or_ maj.rule
    
    ## version.unf defined per default
    ## Tetrad CPC works with version.unf=c(2,1)
    ## see comment on pc.cons.intern for description of version.unf
    pc. <- pc.cons.intern(skel, suffStat, indepTest, alpha,
                          version.unf=c(2,1), maj.rule=maj.rule, verbose=verbose)
    udag2pdagRelaxed(pc.$sk, verbose=verbose,
                     unfVect=pc.$unfTripl, solve.confl=solve.confl)
  }
}
CausaleffectsDirectionsWeijia<-function(miRNA,mRNA,stdData,alp, nmiRs,fixedEdge=NULL,fixedDirections=NULL, solve.confl=FALSE,skel.method=skel.method){
  
  #Learning PCDAG
  suffStat<-list(C=cor(stdData), n=nrow(stdData))  
  
  print('Begin PC')
  pc.fit<-modifiedpcDirectionsWeijia(suffStat, indepTest=gaussCItest, p=ncol(stdData), alpha=alp,fixedEdge=fixedEdge,fixedDirections=fixedDirections,skel.method=skel.method,solve.confl=solve.confl)               
  print("PC finished")
  
  nmRs<-ncol(stdData)-nmiRs
  result<-matrix(nrow=nmRs, ncol=nmiRs)
  #for (l in 1:nmiRs){
  l=1
  while(l<=nmiRs){
    #Inferring causal effects
    caef<-try(idaFast(l,(nmiRs+1):ncol(stdData),cov(stdData), pc.fit@graph ))
    print(l)
    if(class(caef) == "try-error") {result=NA;print("singular error");break}
    else{
      #min of absolute values.
      caef1<-matrix(nrow=nmRs,ncol=1)
      for (k in 1:nmRs){
        caefabs<-abs(caef)
        index<-which(caefabs==min(caefabs[k,]), arr.ind=TRUE)
        pos<-index[1,2]
        caef1[k,]<-caef[k,pos]
      }
      result[,l]<-caef1
      l=l+1
    }
  }
  rownames(result)<-mRNA
  colnames(result)<-miRNA
  print('Causaul Finish!')
  return(result)
}

skeletonThuc <- function(suffStat, indepTest, alpha, labels, p,
                      method = c("stable", "original", "stable.fast"), m.max = Inf,
                      fixedGaps = NULL, fixedEdges = NULL,
                      NAdelete = TRUE, verbose = FALSE)
{
 
  cl <- match.call()
  if(!missing(p)) stopifnot(is.numeric(p), length(p <- as.integer(p)) == 1, p >= 2)
  if(missing(labels)) {
    if(missing(p)) stop("need to specify 'labels' or 'p'")
    labels <- as.character(seq_len(p))
  } else { ## use labels ==> p  from it
    stopifnot(is.character(labels))
    if(missing(p)) {
      p <- length(labels)
    } else if(p != length(labels))
      stop("'p' is not needed when 'labels' is specified, and must match length(labels)")
    else
      message("No need to specify 'p', when 'labels' is given")
  }
  seq_p <- seq_len(p)
  method <- match.arg(method)
  ## C++ version still has problems under Windows; will have to check why
  #  if (method == "stable.fast" && .Platform$OS.type == "windows") {
  #    method <- "stable"
  #    warning("Method 'stable.fast' is not available under Windows; using 'stable' instead.")
  #  }
  
  ## G := !fixedGaps, i.e. G[i,j] is true  iff  i--j  will be investigated
  if (is.null(fixedGaps)) {
    G <- matrix(TRUE, nrow = p, ncol = p)
    diag(G) <- FALSE
  }
  else if (!identical(dim(fixedGaps), c(p, p)))
    stop("Dimensions of the dataset and fixedGaps do not agree.")
  else if (!identical(fixedGaps, t(fixedGaps)) )
    stop("fixedGaps must be symmetric")
  else
    G <- !fixedGaps
  
  if (any(is.null(fixedEdges))) { ## MM: could be sparse
    #print("fixedEdges is null, don't know why")
    fixedEdges <- matrix(rep(FALSE, p * p), nrow = p, ncol = p)
  }
  else if (!identical(dim(fixedEdges), c(p, p)))
    stop("Dimensions of the dataset and fixedEdges do not agree.")
  # else if (fixedEdges != t(fixedEdges))
  #  print("fixedEdges must be symmetric. Don't worry about the warning, just keep going (Thuc)")
  #Thuc added
  #else print(fixedEdges)
  
  if (method == "stable.fast") {
    ## Do calculation in C++...
    if (identical(indepTest, gaussCItest))
      indepTestName <- "gauss"
    else
      indepTestName <- "rfun"
    options <- list(
      verbose = as.integer(verbose), 
      m.max = as.integer(ifelse(is.infinite(m.max), p, m.max)),
      NAdelete = NAdelete)
    res <- .Call("estimateSkeleton", G, suffStat, indepTestName, indepTest, alpha, fixedEdges, options);
    G <- res$amat
    # sepset <- res$sepset
    sepset <- lapply(seq_p, function(i) c(
      lapply(res$sepset[[i]], function(v) if(identical(v, as.integer(-1))) NULL else v),
      vector("list", p - length(res$sepset[[i]])))) # TODO change convention: make sepset triangular
    pMax <- res$pMax
    n.edgetests <- res$n.edgetests
    ord <- length(n.edgetests) - 1
  }
  else {
    ## Original R version
    
    pval <- NULL
    sepset <- lapply(seq_p, function(.) vector("list",p))# a list of lists [p x p]
    ## save maximal p value
    pMax <- matrix(-Inf, nrow = p, ncol = p)
    diag(pMax) <- 1
    done <- FALSE
    ord <- 0
    n.edgetests <- numeric(1)# final length = max { ord}
    while (!done && any(G) && ord <= m.max) {
      n.edgetests[ord1 <- ord+1L] <- 0
      done <- TRUE
      ind <- which(G, arr.ind = TRUE)
      ## For comparison with C++ sort according to first row
      ind <- ind[order(ind[, 1]), ]
      remainingEdgeTests <- nrow(ind)
      cat("Order=", ord, "; remaining edges:", remainingEdgeTests,"\n",sep="")
      if(method == "stable") {
        ## Order-independent version: Compute the adjacency sets for any vertex
        ## Then don't update when edges are deleted
        G.l <- split(G, gl(p,p))
      }
      for (i in 1:remainingEdgeTests) {
        if (verbose && i%%100 == 0) cat("|i=", i, "|iMax=", nrow(ind), "\n")
        x <- ind[i, 1]
        y <- ind[i, 2]
        if (G[y, x] && !fixedEdges[y, x]) {
          #print(G)
          #print(fixedEdges)
          nbrsBool <- if(method == "stable") G.l[[x]] else G[,x]
          nbrsBool[y] <- FALSE
          nbrs <- seq_p[nbrsBool]
          length_nbrs <- length(nbrs)
          if (length_nbrs >= ord) {
            if (length_nbrs > ord)
              done <- FALSE
            S <- seq_len(ord)
            repeat { ## condition w.r.to all  nbrs[S] of size 'ord'
              n.edgetests[ord1] <- n.edgetests[ord1] + 1
              pval <- indepTest(x, y, nbrs[S], suffStat)
              if (verbose)
                cat("x=", x, " y=", y, " S=", nbrs[S], ": pval =", pval, "\n")
              if(is.na(pval))
                pval <- as.numeric(NAdelete) ## = if(NAdelete) 1 else 0
              if (pMax[x, y] < pval)
                pMax[x, y] <- pval
              if(pval >= alpha) { # independent
                G[x, y] <- G[y, x] <- FALSE
                sepset[[x]][[y]] <- nbrs[S]
                break
              }
              else {
                nextSet <- getNextSet(length_nbrs, ord, S)
                if (nextSet$wasLast)
                  break
                S <- nextSet$nextSet
              }
            } ## repeat
          }
        }
      }# for( i )
      ord <- ord + 1
    } ## while()
    for (i in 1:(p - 1)) {
      for (j in 2:p)
        pMax[i, j] <- pMax[j, i] <- max(pMax[i, j], pMax[j,i])
    }
  }
  
  ## transform matrix to graph object :
  Gobject <-
    if (sum(G) == 0) {
      new("graphNEL", nodes = labels)
    } else {
      colnames(G) <- rownames(G) <- labels
      as(G+fixedEdges,"graphNEL")
    }
  
  ## final object
  new("pcAlgo", graph = Gobject, call = cl, n = integer(0),
      max.ord = as.integer(ord - 1), n.edgetests = n.edgetests,
      sepset = sepset, pMax = pMax, zMin = matrix(NA, 1, 1))
  
}## end{ skeleton }




modifiedpcThuc <- function(suffStat, indepTest, alpha, labels, p,
                       fixedGaps = NULL, fixedEdges = NULL, NAdelete = TRUE, m.max = Inf,
                       u2pd = c("relaxed", "rand", "retry"),
                       skel.method = c("stable", "original", "stable.fast"),
                       conservative = FALSE, maj.rule = FALSE,
                       solve.confl = FALSE, verbose = FALSE){
  ## Purpose: Perform PC-Algorithm, i.e., estimate skeleton of DAG given data
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - dm: Data matrix (rows: samples, cols: nodes)
  ## - C: correlation matrix (only for continuous)
  ## - n: sample size
  ## - alpha: Significance level of individual partial correlation tests
  ## - corMethod: "standard" or "Qn" for standard or robust correlation
  ##              estimation
  ## - G: the adjacency matrix of the graph from which the algorithm
  ##      should start (logical)
  ## - datatype: distinguish between discrete and continuous data
  ## - NAdelete: delete edge if pval=NA (for discrete data)
  ## - m.max: maximal size of conditioning set
  ## - u2pd: Function for converting udag to pdag
  ##   "rand": udag2pdag
  ##   "relaxed": udag2pdagRelaxed
  ##   "retry": udag2pdagSpecial
  ## - gTrue: Graph suffStatect of true DAG
  ## - conservative: If TRUE, conservative PC is done
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 26 Jan 2006; Martin Maechler
  ## Modifications: Sarah Gerster, Diego Colombo, Markus Kalisch
  
  ## Initial Checks
  cl <- match.call()
  if(!missing(p)) stopifnot(is.numeric(p), length(p <- as.integer(p)) == 1, p >= 2)
  if(missing(labels)) {
    if(missing(p)) stop("need to specify 'labels' or 'p'")
    labels <- as.character(seq_len(p))
  } else { ## use labels ==> p  from it
    stopifnot(is.character(labels))
    if(missing(p)) {
      p <- length(labels)
    } else if(p != length(labels))
      stop("'p' is not needed when 'labels' is specified, and must match length(labels)")
    else
      message("No need to specify 'p', when 'labels' is given")
  }
  seq_p <- seq_len(p)
  
  u2pd <- match.arg(u2pd)
  skel.method <- match.arg(skel.method)
  if(u2pd != "relaxed") {
    if (conservative || maj.rule)
      stop("Conservative PC and majority rule PC can only be run with 'u2pd = relaxed'")
    
    if (solve.confl)
      stop("Versions of PC using lists for the orientation rules (and possibly bi-directed edges)\n can only be run with 'u2pd = relaxed'")
  }
  
  if (conservative && maj.rule) stop("Choose either conservative PC or majority rule PC!")
  
  ## Skeleton
  skel <- skeletonThuc(suffStat, indepTest, alpha, labels=labels, method = skel.method,
                    fixedGaps=fixedGaps, fixedEdges=fixedEdges,
                    NAdelete=NAdelete, m.max=m.max, verbose=verbose)
  skel@call <- cl # so that makes it into result
  
  ## Orient edges
  if (!conservative && !maj.rule) {
    switch (u2pd,
            "rand" = udag2pdag(skel),
            "retry" = udag2pdagSpecial(skel)$pcObj,
            "relaxed" = udag2pdagRelaxed(skel, verbose=verbose, solve.confl=solve.confl))
  }
  else { ## u2pd "relaxed" : conservative _or_ maj.rule
    
    ## version.unf defined per default
    ## Tetrad CPC works with version.unf=c(2,1)
    ## see comment on pc.cons.intern for description of version.unf
    pc. <- pc.cons.intern(skel, suffStat, indepTest, alpha,
                          version.unf=c(2,1), maj.rule=maj.rule, verbose=verbose)
    udag2pdagRelaxed(pc.$sk, verbose=verbose,
                     unfVect=pc.$unfTripl, solve.confl=solve.confl)
  }
}

lm.cov <- function (C, y, x) {
  solve(C[x, x], C[x, y, drop = FALSE])[1, ]
}

has.new.coll <- function(amat,amatSkel, x, pa1, pa2.t, pa2.f) {
  ## Check if undirected edges that are pointed to x create a new v-structure
  ## Additionally check, if edges that are pointed away from x create
  ## new v-structure; i.e. x -> pa <- papa would be problematic
  ## pa1 are definit parents of x
  ## pa2 are undirected "parents" of x
  ## pa2.t are the nodes in pa2 that are directed towards pa2
  ## pa2.f are the nodes in pa2 that are directed away from pa2
  ## Value is TRUE, if new collider is introduced
  res <- FALSE
  if (length(pa2.t) > 0 && !all(is.na(pa2.t))) {
    ## check whether all pa1 and all pa2.t are connected;
    ## if no, there is a new collider
    if (length(pa1) > 0 && !all(is.na(pa1))) {
      res <- min(amatSkel[pa1, pa2.t]) == 0 ## TRUE if new collider
    }
    ## in addition, all pa2.t have to be connected
    if (!res && length(pa2.t) > 1) {
      A2 <- amatSkel[pa2.t,pa2.t]
      diag(A2) <- 1
      res <- min(A2) == 0 ## TRUE if new collider
    }
  }
  if (!res && length(pa2.f) > 0 && !all(is.na(pa2.f))) {
    ## consider here only the DIRECTED Parents of pa2.f
    ## remove undirected edges
    A <- amat-t(amat)
    A[A<0] <- 0
    ## find parents of pa2.f
    cA <- colSums(A[pa2.f,,drop=FALSE])
    papa <- setdiff(which(cA != 0), x)
    ## if any node in papa is not directly connected to x, there is a new
    ## collider
    if (length(papa) > 0)
      res <- min(amatSkel[x,papa]) == 0 ## TRUE if new collider
  }
  res
}
wgtMatrix.0 <- function(g, transpose = TRUE){
  ## Purpose: work around "graph" package's  as(g, "matrix") bug
  ## ----------------------------------------------------------------------
  ## ACHTUNG: mat_[i,j]==1 iff j->i,
  ## whereas with as(g,"matrix") mat_[i,j]==1 iff i->j
  ## ----------------------------------------------------------------------
  ## Arguments: g: an object inheriting from (S4) class "graph"
  ## ----------------------------------------------------------------------
  ## Author: Martin Maechler, based on Seth Falcon's code;  Date: 12 May 2006
  
  ## MM: another buglet for the case of  "no edges":
  if(numEdges(g) == 0) {
    p <- length(nd <- nodes(g))
    return( matrix(0, p,p, dimnames = list(nd, nd)) )
  }
  ## Usual case, when there are edges:
  if(!("weight" %in% names(edgeDataDefaults(g))))
    edgeDataDefaults(g, "weight") <- 1L
  w <- unlist(edgeData(g, attr = "weight"))
  ## we need the *transposed* matrix typically:
  tm <- if(transpose) t(as(g, "matrix")) else as(g, "matrix")
  ## now is a 0/1 - matrix (instead of 0/wgts) with the 'graph' bug
  if(any(w != 1)) ## fix it
    tm[tm != 0] <- w
  ## tm_[i,j]==1 iff i->j
  tm
}
