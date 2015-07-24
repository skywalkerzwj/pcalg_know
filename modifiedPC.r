library(pcalg)
mskeleton <- function(suffStat, indepTest, alpha, labels, p,
		     method = c("stable", "original", "stable.fast"), m.max = Inf,
		     fixedGaps = NULL, fixedEdges = NULL,
		     NAdelete = TRUE, verbose = FALSE)
{ cl <- match.call()
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
	#print("fixedEdges is null, don't know why")
    fixedEdges <- matrix(rep(FALSE, p * p), nrow = p, ncol = p)
  }
  else if (!identical(dim(fixedEdges), c(p, p)))
    stop("Dimensions of the dataset and fixedEdges do not agree.")
 # else if (fixedEdges != t(fixedEdges))
  #  print("fixedEdges must be symmetric. Don't worry about the warning, just keep going (Thuc)")
#Thuc added
  #else print(fixedEdges)
	

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
      if (verbose)
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

skeleton <- function(suffStat, indepTest, alpha, labels, p,
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
 
  G <- matrix(TRUE, nrow = p, ncol = p)
  diag(G) <- FALSE
 
  if (any(is.null(fixedEdges))) { ## MM: could be sparse
    fixedEdges <- matrix(rep(FALSE, p * p), nrow = p, ncol = p)
  }
  else if (!identical(dim(fixedEdges), c(p, p)))
    stop("Dimensions of the dataset and fixedEdges do not agree.")
#   else if (fixedEdges != t(fixedEdges))
#     stop("fixedEdges must be symmetric")
  
    ## Original R version    
    pval <- NULL
    ## save maximal p value
    done <- FALSE
    ord <- 0
    while (!done && any(G) && ord <= m.max) { #any(G-fixedEdge)
      done <- TRUE      
#       test=G-fixedEdges
#       test=test!=0
#       ind <- which(test, arr.ind = TRUE)  #modified 16/10
      ind <- which(G, arr.ind = TRUE)
      ## For comparison with C++ sort according to first row
      ind <- ind[order(ind[, 1]), ]
      remainingEdgeTests <- nrow(ind)
        cat("Order=", ord, "; remaining edges:", remainingEdgeTests,"\n",sep="")
      
      for (i in 1:remainingEdgeTests) {
#         if (ord>=4) print(i)
        x <- ind[i, 1]
        y <- ind[i, 2]
        if (G[y, x] && !fixedEdges[y, x]) {
          nbrsBool <- G[,x]
          nbrsBool[y] <- FALSE
          nbrs <- seq_p[nbrsBool]
          length_nbrs <- length(nbrs)
          if (length_nbrs >= ord) {
            if (length_nbrs > ord)
              done <- FALSE
            S <- seq_len(ord)
            repeat { ## condition w.r.to all  nbrs[S] of size 'ord'
              pval <- indepTest(x, y, nbrs[S], suffStat)
              if(is.na(pval))
                pval <- as.numeric(NAdelete) ## = if(NAdelete) 1 else 0
              if(pval >= alpha) { # independent
#                 G[x, y] <- G[y, x] <- FALSE
                G[y, x] <- FALSE
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

  ## transform matrix to graph object :
  Gobject <-
    if (sum(G) == 0) {
      new("graphNEL", nodes = labels)
    } else {
      colnames(G) <- rownames(G) <- labels
      as(G,"graphNEL")
    }
  
  ## final object
  new("pcAlgo", graph = Gobject, call = cl, n = integer(0),
      max.ord = as.integer(ord - 1), zMin = matrix(NA, 1, 1))
#   new("pcAlgo", graph = Gobject, call = cl, n = integer(0),
#     max.ord = as.integer(ord - 1), n.edgetests = n.edgetests,
#     sepset = sepset, pMax = pMax, zMin = matrix(NA, 1, 1))
#   
}## end{ skeleton }

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
    beta.tmp <- lm.cov(mcov,y.pos.set,c(x.pos,pa1)) ####
    beta.tmp[y.pos.set %in% pa1] <- 0
    beta.hat <- cbind(beta.tmp)
  } else {    ## at least one undirected parent
    ## no member of pa2
    pa2.f <- pa2
    pa2.t <- NA
    beta.hat <-
      if (!has.new.coll(amat,amatSkel,x.pos,pa1,pa2.t,pa2.f)) {
        beta.tmp <- lm.cov(mcov,y.pos.set,c(x.pos,pa1)) ####
        beta.tmp[y.pos.set %in% pa1] <- 0
        cbind(beta.tmp)
      } # else NULL
    
    ## exactly one member of pa2
    for (i2 in seq_along(pa2)) {
      pa2.f <- pa2[-i2]
      pa2.t <- pa2[i2]
      if (!has.new.coll(amat,amatSkel,x.pos,pa1,pa2.t,pa2.f)) {
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


wgtMatrix.0 <- function(g, transpose = TRUE)
{
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

modifiedpc <- function(suffStat, indepTest, alpha, labels, p,
               fixedGaps = NULL, fixedEdges = NULL, NAdelete = TRUE, m.max = Inf,
               u2pd = c("relaxed", "rand", "retry"),
               skel.method = c("stable", "original", "stable.fast"),
               conservative = FALSE, maj.rule = FALSE,
               solve.confl = FALSE, verbose = FALSE)
{
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
  skel <- skeleton(suffStat, indepTest, alpha, labels=labels, method = skel.method,
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
