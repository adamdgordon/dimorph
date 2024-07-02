#' Get the Structure Component for Dimorphism Resampling
#' 
#' Given a list of metric datasets, calculates the shared data structure that allows
#' for all samples to be included in a common resampling analysis.  This generates
#' the \code{struc} component to be included in various functions in the \code{dimorph}
#' package.
#' @seealso \code{\link{SSDtest}}
#' @examples
#' SSDvars <- c("FHSI", "TPML", "TPMAP", "TPLAP", "HHMaj",
#'              "HHMin", "RHMaj", "RHMin", "RDAP", "RDML")
#' getstructure(list(apelimbart[apelimbart$Species=="Gorilla gorilla", SSDvars],
#'                   apelimbart[apelimbart$Species=="Homo sapiens", SSDvars]))
#' getstructure(list(apelimbart[apelimbart$Species=="Gorilla gorilla", SSDvars],
#'                   apelimbart[apelimbart$Species=="Homo sapiens", SSDvars],
#'                   fauxil[fauxil$Species=="Fauxil sp. 1", SSDvars]))
#' getstructure(list(apelimbart[apelimbart$Species=="Gorilla gorilla", SSDvars],
#'                   apelimbart[apelimbart$Species=="Homo sapiens", SSDvars],
#'                   fauxil[fauxil$Species=="Fauxil sp. 1", SSDvars],
#'                   fauxil[fauxil$Species=="Fauxil sp. 2", SSDvars]))
#' getstructure(list(fauxil[fauxil$Species=="Fauxil sp. 1", SSDvars]))
#' getstructure(list(fauxil[fauxil$Species=="Fauxil sp. 2", SSDvars]))
#' getstructure(list(fauxil[fauxil$Species=="Fauxil sp. 2", SSDvars]), forcematrix=T)
#' @export
getstructure <- function(x, forcematrix=F) {
  if(class(x)!="list") stop("'x' must be a list of dataframes or matrices of samples to generate a shared structure from.")
  ndat <- length(x)
  if (ndat < 2) {
	struc <- x[[1]]
	struc <- struc/struc
	nspecpervar <- apply(struc, 2, function(x) sum(!is.na(x)))
	goodvars <- nspecpervar > 1
	struc <- struc[,goodvars, drop=F]
	nspecpervar <- nspecpervar[goodvars]
	# forcematrix prevents matrix structure from being converted to vector structure
	if (!forcematrix) if (prod(nspecpervar==nspecpervar[1])==0) struc <- nspecpervar
	#stop("'getstructure' is only needed when comparing two or more datasets.")
  } # end ndat==1
  else { # if ndat > 1
  gooddat <- unlist(lapply(x, function(x) class(x)=="data.frame" | class(x)=="matrix"))
  if (prod(gooddat)==0) stop("Each element of 'x' must be a dataframe or matrix.")
  # first match variables
  sharedvars <- colnames(x[[1]])
  for (i in 2:ndat) {
    sharedvars <- sharedvars[sharedvars %in% colnames(x[[i]])]
  }
  if (length(sharedvars)==0) stop("The datasets in 'x' must share at least one variable in common.")
  for (i in 1:ndat) x[[i]] <- x[[i]][,sharedvars,drop=F]
  # remove any rows missing all data
  for (i in 1:ndat) x[[i]] <- x[[i]][apply(x[[i]], 1, function (x) return(sum(is.na(x)))) < length(sharedvars),,drop=F]
  # now identify complete and incomplete datasets
  complete <- unlist(lapply(x, function(x) as.logical(prod(complete.cases(x)))))
  if (prod(complete)==0) { # if there are some datasets with missing data
    if (sum(!complete) > 1) { # more than one dataset with missing data
	  nspecpervar <- do.call(rbind, lapply(x, function(x) apply(x, 2, function(x) sum(!is.na(x)))))
	  nspecpervar <- apply(nspecpervar, 2, min)
	  goodvars <- nspecpervar > 1
	  if (sum(goodvars)==0) stop("At least one variable must have at least two measurements in all datasets.")
	  if (prod(goodvars)==0) {
	    # remove variables with fewer than two measurements in x
	    sharedvars <- sharedvars[goodvars]
		for (i in 1:ndat) x[[i]] <- x[[i]][,sharedvars,drop=F]
		# go back through x removing any rows that now have no data
	    for (i in 1:ndat) x[[i]] <- x[[i]][apply(x[[i]], 1, function (x) return(sum(is.na(x)))) < length(sharedvars),,drop=F]
	    # recalculate nspecpervar
		nspecpervar <- do.call(rbind, lapply(x, function(x) apply(x, 2, function(x) sum(!is.na(x)))))
		nspecpervar <- apply(nspecpervar, 2, min)
	  }
	  struc <- nspecpervar
	}
	else { # only one dataset with missing data
	  struc <- x[[which(!complete)]]
	  struc <- struc/struc
	  nspecpervar <- apply(struc, 2, function(x) sum(!is.na(x)))
	  goodvars <- nspecpervar > 1 # need at least two measurements in each included variable
	  if (sum(goodvars)==0) stop("At least one variable must have at least two measurements in all datasets.")
	  if (prod(goodvars)==0) {
	    # remove variables with fewer than two measurements in x and struc
	    sharedvars <- sharedvars[goodvars]
		struc <- struc[,sharedvars]
		for (i in 1:ndat) x[[i]] <- x[[i]][,sharedvars,drop=F]
		# go back through x and struc removing any rows that now have no data
	    for (i in 1:ndat) x[[i]] <- x[[i]][apply(x[[i]], 1, function (x) return(sum(is.na(x)))) < length(sharedvars),,drop=F]
	    struc <- struc[apply(struc, 1, function (x) return(sum(is.na(x)))) < length(sharedvars),,drop=F]
	  }
	  nspecstruc <- nrow(struc)
	  nspecadequate <- unlist(lapply(x, function(x) nrow(x) >= nspecstruc))
	  if (prod(nspecadequate)==0) stop("Not all of the complete datasets have large enough sample sizes to sample down to the missing data structure.")
	}
  } # end missing data section
  else { # if all datasets are complete
    nspec <- unlist(lapply(x, function(x) nrow(x)))
	if (prod(nspec > 1)==0)  stop("Not all of the datasets have at least two specimens.")
	struc <- x[[which(nspec==min(nspec))[1]]]
	struc <- struc/struc
  }
  } # end ndat > 1
  return(struc)
}


#' Get Resampling Addresses for Comparative Sample for Estimating Dimorphism
#' 
#' Function called by other functions in package \code{'dimorph'} to generate a set of resampled addresses in a 
#'   dataset for use in resampling analyses of dimorphism.
#' @param comparative A matrix or data frame of measurements from a comparative sample, with rows
#'   corresponding to individual specimens and columns corresponding to size variables.  Sex data
#'   should not be included.  If \code{comparative} contains \code{NA}s then \code{struc} must be supplied
#'   as a vector rather than as a matrix or data frame.
#' @param struc Structure information for the data set being compared against, typically one or 
#'   more fossil samples with missing data. \code{struc} must be either (1) a matrix or dataframe
#'   of measurements (which can include \code{NA}s) or (2) a vector of integer sample sizes for each
#'   variable. 
#' @param compsex A vector indicating sex for the individuals in \code{comparative}.  Sex
#'   information is not included in any calculations in this function but will be included
#'   as a list element in the returned object.  Defaults to \code{NULL}.
#' @param exact Logical scalar specifying whether or not to sample exactly once all possible 
#'   unique combinations of the comparative data matching the pattern in \code{struc}.  This 
#'   procedure takes into account any missing data pattern in \code{struc}.  For example, for 
#'   three observations, "ABC" and "CAB" are identical samples if there is no missing data, but
#'   they are different samples if the first observation has a missing data pattern that differs
#'   from that of the second and third observations. \code{exact} defaults 
#'   to \code{FALSE}.  If \code{FALSE}, Monte Carlo sampling is used instead.
#' @param limit An upper limit for the number of unique samples allowable for exact resampling.
#'   if \code{exact} is \code{TRUE} and the number of unique combinations would exceed \code{limit} 
#'   then Monte Carlo sampling is used instead.  \code{limit} is ignored if \code{exact} is \code{FALSE}.
#' @param nResamp Integer specifying the number of resampling iterations to use when \code{exact} is
#'   \code{FALSE} or when \code{exact} is \code{TRUE} but the total number of unique combinations exceeds
#'   \code{limit}.
#' @param matchvars Logical scalar specifying whether to compare the variable names in
#'   \code{comparative} and \code{struc} and pare them both down to the set of shared variables. 
#'   If \code{FALSE} and variable names differ then an error will be returned.  Defaults to
#'   \code{FALSE}.
#' @param replace Logical scalar passed to \code{\link[base]{sample}} specifying whether or not
#'   to sample with replacement.  Defaults to \code{FALSE}.
#' @return A list of class \code{dimorphAds} containing resampled addresses and the information used to 
#'   generate them.
#' @seealso \code{\link{resampleSSD}}
#' @examples
#' ## Univariate addresses
#' SSDvars <- c("HHMaj")
#' addressesUni <- getsampleaddresses(comparative=apelimbart[apelimbart$Species=="Gorilla gorilla", SSDvars, drop=F],
#'                                    struc=fauxil[fauxil$Species=="Fauxil sp. 1", SSDvars, drop=F],
#'                                    compsex=apelimbart[apelimbart$Species=="Gorilla gorilla", "Sex"],
#'                                    exact=T, matchvars=T, replace=F)
#' addressesUni
#' str(addressesUni)
#' 
#' ## Multivariate addresses
#' SSDvars <- c("HHMaj","RHMaj","FHSI","TPML")
#' addressesMulti1 <- getsampleaddresses(comparative=apelimbart[apelimbart$Species=="Gorilla gorilla", SSDvars],
#'                                       struc=fauxil[fauxil$Species=="Fauxil sp. 1", SSDvars],
#'                                       compsex=apelimbart[apelimbart$Species=="Gorilla gorilla", "Sex"],
#'                                       nResamp=10000, matchvars=T, replace=F)
#' addressesMulti1
#' str(addressesMulti1)
#' 
#' addressesMulti2 <- getsampleaddresses(comparative=apelimbart[apelimbart$Species=="Gorilla gorilla", SSDvars],
#'                                       struc=setNames(c(2,2,4,5), SSDvars),
#'                                       compsex=apelimbart[apelimbart$Species=="Gorilla gorilla", "Sex"],
#'                                       nResamp=10000, matchvars=T, replace=F)
#' addressesMulti2
#' str(addressesMulti2)
#' 
#' addressesMulti3 <- getsampleaddresses(comparative=apelimbart[apelimbart$Species=="Gorilla gorilla", "FHSI"],
#'                                       struc=setNames(4, "FHSI"),
#'                                       compsex=apelimbart[apelimbart$Species=="Gorilla gorilla", "Sex"],
#'                                       nResamp=10000, matchvars=T, replace=F)
#' addressesMulti3
#' str(addressesMulti3)
#' 
#' # Now with missing data also in the comparative sample
#' SSDvars <- c("FHSI", "TPML", "TPMAP", "TPLAP", "HHMaj",
#'              "HHMin", "RHMaj", "RDAP", "RDML")
#' Fs1 <- fauxil[fauxil$Species=="Fauxil sp. 1", SSDvars]
#' Fs2 <- fauxil[fauxil$Species=="Fauxil sp. 2", SSDvars]
#' n.Fs1 <- apply(Fs1, 2, function(x) sum(!is.na(x)))
#' n.Fs2 <- apply(Fs2, 2, function(x) sum(!is.na(x)))
#' n.Fs1
#' n.Fs2
#' n.min <- apply(rbind(n.Fs1,n.Fs2), 2, min)
#' gorAds <- getsampleaddresses(comparative=apelimbart[apelimbart$Species=="Gorilla gorilla", SSDvars],
#'                              struc=Fs1,
#'                              exact=T, limit=500000, replace=F)
#' gorAds
#' str(gorAds)
#' 
#' Fs1Ads <- getsampleaddresses(comparative=Fs1, struc=n.min,
#'                              exact=T, limit=500000, replace=F)
#' Fs1Ads
#' str(Fs1Ads)
#' 
#' Fs1AdsReplace <- getsampleaddresses(comparative=Fs1, struc=n.min,
#'                                     exact=T, limit=500000, nResamp=1000,
#'                                     replace=T)
#' Fs1AdsReplace
#' str(Fs1AdsReplace)
#' 
#' Fs2Ads <- getsampleaddresses(comparative=Fs2, struc=n.min,
#'                              exact=T, limit=500000, replace=F)
#' Fs2Ads
#' str(Fs2Ads)
#' @export
getsampleaddresses <- function(comparative,
                               struc,
							   compsex=NULL,
							   sex.female=1,
							   nResamp=10000,
							   exact=F,
							   limit=500000,
                               matchvars=F,
							   replace=F) {
  # "struc" is either a matrix of measurement data (can have missing data) or a 
  # vector corresponding to sample sizes for each variable
  if (class(comparative)=="numeric" | class(comparative)=="integer") comparative <- data.frame(VAR=comparative)
  { # if sex is present, make sure it matches nrow of comparative and set flag according to forcemixed
    if(!is.null(compsex)) {
      if (!length(compsex)==nrow(comparative)) stop("If present, 'compsex' must have a number of elements equal to the number of individuals in 'comparative'.")
      compsex <- droplevels(as.factor(compsex))
      if (!(nlevels(compsex)==1|nlevels(compsex)==2)) stop("The sex variable can have no more than two levels.")
      if (!(sex.female==1 | sex.female==2)) stop("The level in 'compsex' which corresponds to females must be indicated.")
    }
	else sex.female <- NULL
  }
  if (!(is.data.frame(comparative)|is.matrix(comparative))) stop("Argument 'comparative' must be a vector, data frame, or matrix.")
#  if(prod(complete.cases(comparative))==0) warning("The comparative data set is not complete for all specimens.  Resampling addresses will not be affected, but these addresses may generate errors when used by 'resampleSSD()'.")
  if(prod(complete.cases(comparative))==0 & !(is.vector(struc))) stop("The comparative data set is not complete for all specimens.  'struc' should be provided as a vector of sample sizes per variable rather than a data matrix.  Type 'help(getsampleaddresses)' to see example.")
  varsComp <- colnames(comparative)
  {
    if (is.vector(struc)) {
      snames <- names(struc)
      struc <- as.integer(struc)
      names(struc) <- snames
      if(!class(struc)=="integer") stop("The structure object must either be (1) a matrix or dataframe of measurements (which can include NAs) or (2) a vector of integer sample sizes for each variable.")
      varsStruc <- names(struc)
	  if (is.null(varsStruc)) {
	    {if (length(struc)==1) names(struc) <- "VAR"
		else names(struc) <- paste0("VAR_", 1:length(struc))}
		varsStruc <- names(struc)
	  }
	  if (length(varsStruc)==1 & length(varsComp)==1) {
	    varsComp <- varsStruc
		colnames(comparative) <- varsComp
	  }
    }
    else if (is.data.frame(struc)|is.matrix(struc)) {
      struc <- data.frame(struc)
	  # delete specimens with no measurements
	  struc <- struc[apply(struc, 1, function(x) !sum(is.na(x))==length(x)), , drop=F] 
      varsStruc <- colnames(struc)
	  if (length(varsStruc)==1 & length(varsComp)==1 & varsComp[1]=="VAR") {
	    varsComp <- varsStruc
		colnames(comparative) <- varsComp
	  }
    }
    else stop("The structure object must either be (1) a matrix or dataframe of measurements (which can include NAs) or (2) a vector of integer sample sizes for each variable.")
  }
  if (matchvars) {
    vars <- varsComp[varsComp %in% varsStruc]
    comparative <- comparative[,vars, drop=FALSE]
    struc <- struc[vars]
    varsComp <- colnames(comparative)
    {
      if (is.vector(struc)) varsStruc <- names(struc)
      else varsStruc <- colnames(struc)
    }
	rm(vars)
  }
  if(!identical(sort(varsComp), sort(varsStruc))) stop("The comparative data set and structure object must have the same variables.")
  out <- NULL
  { # if loop for struc classes
    if (is.vector(struc)) { # vector struc procedure
	  if (length(varsComp)!=length(varsStruc)) stop("The vector 'struc' must have the same number of elements as variables in 'comparative'.")
	  #compStruc <- nrow(comparative) - apply(comparative, 2, sumNA)
	  compStruc <- nrow(comparative) - apply(comparative, 2, function (x) return(sum(is.na(x))))
	  over <- compStruc-struc
	  if (sum(over < 0) > 0) stop("The dataset to sample from does not contain enough observations in each variable.")
	  whichisna <- function(x) which(!is.na(x))
	  goodAds <- apply(comparative, 2, whichisna, simplify=F)
      ads <- NULL
      adlist <- NULL  
      if (exact) {
	    {if (replace) exactvec <- choose(compStruc+struc-1, struc) # combination WITH replacement (n+k-1, choose k)
	    else exactvec <- choose(compStruc, struc)}
        nExact <- prod(exactvec)
        { # calculates exact addresses if nExact is less than limit, otherwise pushes to Monte Carlo
          if (nExact > limit) {
            warning(paste0("The number of possible combinations (", nExact, ") exceeds the user-specified limit. Monte Carlo sampling will be used."))
            exact <- F
          }
          else {
            nResamp <- nExact
            nRvars <- sum(exactvec > 1)
            possibleAds <- goodAds
		    {if (replace) {
		      for (i in 1:length(possibleAds)) possibleAds[[i]] <- t(RcppAlgos::comboGeneral(goodAds[[i]], struc[i], repetition=T))
		    }
		    else {
              for (i in 1:length(possibleAds)) possibleAds[[i]] <- combn(goodAds[[i]], struc[i])
		    }} # end if else
            expand <- function(x) return(1:x)
		    metaAds <- expand.grid(apply(data.frame(x=exactvec), 1, expand))
		    adlist <- goodAds
		    for (i in 1:length(adlist)) adlist[[i]] <- possibleAds[[i]][,metaAds[,i]]
          }
        }
      } # end exact==T section
      if (!exact) { # don't use else if because exact can be reset to FALSE in earlier section to trigger this
        adlist <- goodAds
        for(i in 1:length(adlist)) {
          adlist[[i]] <- apply(data.frame(x=rep(length(adlist[[i]]), nResamp)), 1,
		                       sample, size=struc[i], replace=replace)
          adlist[[i]] <- matrix(goodAds[[i]][adlist[[i]]], nrow(adlist[[i]]), nResamp)
        }
      } # end exact==F section
      {if (length(varsComp)==1) {
	    ads <- adlist[[1]]
		adlist <- NULL
		struc <- data.frame(x=rep(1, struc))
		colnames(struc) <- varsComp
		strucClass <- "data.frame"
	  }
	  else {
        nmax <- sum(struc)
		strucClass <- "vector"
        getlistreplaceF <- function(i) {
          tmp <- NULL
          for (j in 1:length(adlist)) tmp <- c(tmp, adlist[[j]][,i])
          tmp <- unique(tmp)
          ret <- rep(NA, nmax)
          ret[1:length(tmp)] <- tmp
          return(ret)
        }
        getlistreplaceT <- function(i) {
          tmp <- NULL
          for (j in 1:length(adlist)) {
            tmp <- c(tmp, table(adlist[[j]][,i]))
          }
          tmp <- sort(tmp, decreasing=T)
          tmp <- tmp[unique(names(tmp))]
          tmp2 <- NULL
          for (k in 1:length(tmp)) tmp2 <- c(tmp2, rep(as.integer(names(tmp)[k]), tmp[k]))
          ret <- rep(NA, nmax)
          ret[1:length(tmp2)] <- tmp2
          return(ret)
        }
        if (!replace) ads <- apply(data.frame(x=1:nResamp), 1, getlistreplaceF)
        if (replace)  ads <- apply(data.frame(x=1:nResamp), 1, getlistreplaceT)
        #ads <- ads[apply(ads, 1, sumNA) < nResamp,] # remove completely empty rows
        ads <- ads[apply(ads, 1, function (x) return(sum(is.na(x)))) < nResamp,] # remove completely empty rows
	  }}
      out <- list(nResamp=nResamp,
	              addresses=ads,
				  adlist=adlist,
				  comparative=comparative,
				  compsex=compsex,
				  sex.female=sex.female,
				  replace=replace,
				  exact=exact,
				  matchvars=matchvars,
				  struc=struc,
				  strucClass=strucClass)
    } # end vector struc procedure
    else if (is.data.frame(struc)|is.matrix(struc)){ # matrix struc procedure
	  if (length(varsStruc)==1) struc <- struc[complete.cases(struc), , drop=F]
	  {if (exact) {
	    # add exact section here - first count up number of different patterns
	    n <- nrow(comparative)
	    strucpat <- apply(struc/struc, 1, paste0, collapse="|")
	    patorder <- order(strucpat)
	    strucpat2 <- strucpat[patorder]
	    pattable <- table(strucpat2)
	    allsame <- F
	    alldiff <- F
	    npat <- length(pattable)
	    if (npat==1) allsame <- T # if this is true, just get combinations where order doesn't matter
	    if (npat==length(strucpat)) alldiff <- T # if this is true, just get permutations where order does matter
	    #npat
	    nExact <- 1
	    nremain <- n
	    for (i in 1:npat) {
		  k <- as.integer(pattable[i])
		  if (replace) nExact <- nExact * choose(nremain+k-1, k)
		  else nExact <- nExact * choose(nremain, k)
		  nremain <- nremain - k
		}
		#nExact
		{ # calculates exact addresses if nExact is less than limit, otherwise pushes to Monte Carlo
		  if (nExact > limit) {
		    warning(paste0("The number of possible combinations (", nExact, ") exceeds the user-specified limit. Monte Carlo sampling will be used."))
			exact <- F
	    }
		  else {
		    if (allsame) { # if this is true, just get combinations where order doesn't matter
			  {if (replace) ads <- t(RcppAlgos::comboGeneral(nrow(comparative), nrow(struc), repetition=T))
			  else ads <- combn(nrow(comparative), nrow(struc))} # end if else 
			} # end allsame
			else if (alldiff) { # if this is true, just get permutations where order does matter
			  {if (replace) ads <- t(RcppAlgos::permuteGeneral(nrow(comparative), nrow(struc), repetition=T))
			  else ads <- t(RcppAlgos::permuteGeneral(nrow(comparative), nrow(struc), repetition=F))} # end if else
			} # end allsame
			else if (!allsame & !alldiff) {
			  adslist <- list(NULL)
			  nremain <- n
			  for (i in 1:npat) {
			    k <- as.integer(pattable[i])
				{if (replace) adslist[[i]] <- t(RcppAlgos::comboGeneral(nremain, k, repetition=T))
				else adslist[[i]] <- combn(nremain, k)} # end if else
				nremain <- nremain - k
			  }
			  # now need to convert to actual addresses
			  convertads <- function(x, newads) {
			    tmp <- (1:n)[-x]
				new <- matrix(tmp[newads], nrow(newads), ncol(newads), byrow=F)
				old <- matrix(x, length(x), ncol(new), byrow=F)
				return(rbind(old, new))
			  }
			  ads <- adslist[[1]]
			  for (i in 2:length(adslist)) {
			    ads <- do.call(cbind, apply(ads, 2, convertads, newads=adslist[[i]], simplify=F))
			  }
			  rm(adslist)
			  # ads is in order of pattern rankings - need to get back to original order
			  ads <- ads[order(patorder),]
			} # end !allsame & !alldiff
			nResamp <- ncol(ads)
		  }
		} # end if else
	  } # end exact section
	  if (!exact) { # don't use else if because exact can be reset to FALSE in earlier section to trigger this
	    ads <- apply(data.frame(x=rep(nrow(comparative),nResamp)), 1, sample, size=nrow(struc), replace=replace)
	  }} # end exact if else
      comparative <- data.frame(comparative)
      struc <- data.frame(struc)
      struc <- struc/struc
	  strucClass <- "data.frame"
      out <- list(nResamp=nResamp,
	              addresses=ads,
				  adlist=NULL,
				  comparative=comparative,
				  compsex=compsex,
				  sex.female=sex.female,
				  replace=replace,
				  exact=exact,
				  matchvars=matchvars,
				  struc=struc,
				  strucClass=strucClass)
    } # end matrix struc procedure
  } # end if loop for struc classes
  class(out) <- "dimorphAds"
  return(out)
}


#' Resample Univariate or Multivariate Dimorphism
#' 
#' Function to generate a set of resampled dimorphism estimates for a univariate or multivariate sample.
#'   Called by \code{\link{SSDtest}}.
#' @param x A matrix or data frame of measurements from a comparative sample, with rows
#'   corresponding to individual specimens and columns corresponding to size variables.  Sex data
#'   should not be included.  Should not include \code{NA}s; resampling addresses will not be
#'   affected, but these addresses may generate errors when used by other functions.
#' @param struc Structure information for the data set being compared against, typically from one or 
#'   more fossil samples with missing data. \code{struc} must be either (1) a matrix or dataframe
#'   of measurements (which can include \code{NA}s) or (2) a vector of integer sample sizes for each
#'   variable. 
#' @param compsex A vector indicating sex for the individuals in \code{x}.  Sex
#'   information is not included in any calculations in this function but will be included
#'   as a list element in the returned object.  Defaults to \code{NULL}.
#' @param npersample Integer specifying the sample size of resampled datasets.  Defaults to \code{NA}.
#'   If set to \code{NA}, the sample size of resampled datasets is set equal to the sample size of 
#'   \code{x}.
#' @param nResamp Integer specifying the number of resampling iterations to calculate addresses
#'   for if Monte Carlo sampling is used.
#' @param exact Logical scalar specifying whether to sample all unique combinations of sample
#'   size \code{npersample} from \code{x}.  Defaults to \code{FALSE}.  If set to \code{FALSE}, or if set 
#'   to \code{TRUE} and the number of unique combinations exceeds \code{limit}, then Monte Carlo sampling
#'   is used instead.
#' @param limit Integer setting the upper limit on the number of unique combinations allowable 
#'   for exact resampling.  If exact resampling would produce more resampled datasets than this number, 
#'   Monte Carlo resampling is used instead.  Defaults to 500,000.
#' @param matchvars Logical scalar specifying whether to compare the variable names in
#'   \code{comparative} and \code{struc} and pare them both down to the set of shared variables. 
#'   If \code{FALSE}  and variable names differ then an error will be returned.  Defaults to
#'   \code{FALSE}.
#' @param replace Logical scalar passed to \code{\link[base]{sample}} specifying whether or not
#'   to sample with replacement.  Defaults to \code{FALSE}.
#' @param datastruc If multivariate data are used, this is a character string specifiying whether to 
#'   incorporate the missing data structure
#'   into dimorphism estimates (\code{"missing"}), whether to downsample to the missing data sample size but 
#'   keep all metric data for the comparative sample (\code{"complete"}), or to perform both types of
#'   resampling separately (\code{"both"}).  Defaults to \code{"missing"}, and ignored if only univariate
#'   data are provided.
#' @param methsMulti A character string specifying the multivariate method used to 
#'   calculate or estimate dimorphism.  Note that regardless of the value of this argument,
#'   multivariate estimation procedures will only be carried out if \code{x} is a multivariate 
#'   dataset.  See \code{\link{dimorph}} for options.
#' @param methsUni A character string specifying the univariate method used to calculate or estimate dimorphism. 
#'   See \code{\link{dimorph}} for options.
#' @param sex.female An integer scalar (1 or 2) specifying which level of \code{sex} 
#'   corresponds to female.  Ignored if \code{sex} is \code{NULL}.  Defaults to 1.
#' @param center A character string specifying the method used to calculate a mean, either \code{"geomean"} 
#'   (default) which uses the geometric mean, or \code{"mean"} which uses the arithmetic mean.  More broadly, 
#'   \code{"geomean"} indicates analyses are conducted in logarithmic data space and \code{"mean"} indicates 
#'   analyses are conducted in raw data space.  Some methods can only be applied in one domain or the other: 
#'   \code{"CV"} and \code{"CVsex"} are always calculated in raw data space and \code{center} will be set to 
#'   \code{"mean"} for these methods regardless of the value set by the user; \code{"MoM"}, \code{"sdlog"}, 
#'   and \code{"sdlogsex"} are always calculated in logarithmic data space and \code{center} will be set to 
#'   \code{"geomean"} for these methods regardless of the value set by the user.
#' @param templatevar A character object or integer value specifying the name or column number 
#'   of the variable in \code{x} to be estimated using the template method.  Ignored if template method 
#'   is not used.  Defaults to \code{NULL}.
#' @param na.rm A logical scalar indicating whether NA values should be stripped before
#'   the computation proceeds.  Defaults to \code{TRUE}.
#' @param ncorrection A logical scalar indicating whether to apply Sokal and Braumann's (1980) 
#'   size correction factor to CV estimates.  Defaults to \code{FALSE}.
#' @param details A logical scalar indicating whether variable name and specimen names 
#'   should be retained (if available) as attributes in the output object.  Defaults to 
#'   \code{FALSE}.
#' @return A list of class \code{dimorphResampledUni} or \code{dimorphResampledMulti} containing a dataframe
#'   with resampled dimorphism estimates and a \code{dimorphAds} object containg resampled addresses produced
#'   by \code{\link{getsampleaddresses}}.  Plotting this object produces violin plots for
#'   all resampled distributions.
#' @seealso \code{\link{bootdimorph}}, \code{\link{dimorph}}, \code{\link{getsampleaddresses}}
#' @examples
#' ## Univariate
#' data(apelimbart)
#' gor <- apelimbart[apelimbart$Species=="Gorilla gorilla",]
#' # this is effectively a bootstrap, although see 'bootdimorph'
#' gorSSD <- resampleSSD(gor[,"FHSI", drop=F], methsUni=c("SSD", "MMR", "BDI"),
#'                       compsex=gor$Sex, nResamp=100, replace=T)
#' gorSSD
#' plot(gorSSD)
#' 
#' # now downsample to fossil sample size and sample without replacement
#' SSDvars <- c("HHMaj")
#' gorSSD1 <- resampleSSD(x=apelimbart[apelimbart$Species=="Gorilla gorilla", SSDvars, drop=F],
#'                        struc=fauxil[fauxil$Species=="Fauxil sp. 1", SSDvars, drop=F],
#'                        compsex=apelimbart[apelimbart$Species=="Gorilla gorilla", "Sex"],
#'                        exact=T, matchvars=T, replace=F, methsUni=c("SSD", "MMR", "BDI"))
#' gorSSD1
#' plot(gorSSD1)
#' 
#' # or run 'getsampleaddresses' first
#' addressesUni <- getsampleaddresses(comparative=apelimbart[apelimbart$Species=="Gorilla gorilla", SSDvars, drop=F],
#'                                    struc=fauxil[fauxil$Species=="Fauxil sp. 1", SSDvars, drop=F],
#'                                    compsex=apelimbart[apelimbart$Species=="Gorilla gorilla", "Sex"],
#'                                    exact=T, matchvars=T, replace=F)
#' gorSSD2 <- resampleSSD(addressesUni, methsUni=c("SSD", "MMR", "BDI"))
#' gorSSD2
#' plot(gorSSD2)
#' 
#' ## Multivariate
#' SSDvars <- c("HHMaj","RHMaj","FHSI","TPML")
#' gorSSDmulti1 <- resampleSSD(x=apelimbart[apelimbart$Species=="Gorilla gorilla", SSDvars],
#'                             struc=fauxil[fauxil$Species=="Fauxil sp. 1", SSDvars],
#'                             compsex=apelimbart[apelimbart$Species=="Gorilla gorilla", "Sex"],
#'                             nResamp=100,
#'                             datastruc="both",
#'                             methsMulti = c("GMM"),
#'                             methsUni = c("SSD", "MMR", "BDI"),
#'                             matchvars=T,
#'                             replace=F)
#' gorSSDmulti1
#' plot(gorSSDmulti1)
#' 
#' # or run 'getsampleaddresses' first
#' addresses <- getsampleaddresses(comparative=
#'                apelimbart[apelimbart$Species=="Gorilla gorilla", SSDvars],
#'                struc=fauxil[fauxil$Species=="Fauxil sp. 1", SSDvars],
#'                compsex=apelimbart[apelimbart$Species=="Gorilla gorilla", "Sex"],
#'                nResamp=100, matchvars=T, replace=F)
#' gorSSDmulti2 <- resampleSSD(x=addresses,
#'                             datastruc="both",
#'                             methsMulti = c("GMM"),
#'                             methsUni = c("SSD", "MMR", "BDI"))
#' gorSSDmulti2
#' plot(gorSSDmulti2)
#' @export
resampleSSD <- function(x,
                        struc=NULL,
						compsex=NULL,
						npersample=NA,
						nResamp=10000,
						exact=F,
						limit=500000,
						matchvars=F,
						replace=F,
						datastruc=NULL, #"missing"
						methsMulti=NULL, # c("GMM")
						methsUni = c("MMR", "BDI"),
						sex.female = 1,
						center = "geomean",
						templatevar = NULL,
						na.rm = T,
						ncorrection=F,
						details = F) {
  if (class(x)=="dimorphAds") {
    if (dim(x$comparative)[2] == 1) { # univariate
	  out <- dimorph::resampleSSDuni(x,
	                                 npersample=npersample,
									 methsUni=methsUni,
									 sex.female=sex.female,
									 center=center,
									 na.rm=na.rm,
									 ncorrection=ncorrection,
									 details=details)
	  out$estimates <- dimorph::logandcalcbiasUni(out$estimates)
	}
	else { #multivariate
	  if (is.null(datastruc)) datastruc <- "missing"
	  if (is.null(methsMulti)) methsMulti <- c("GMM")
	  out <- dimorph::resampleSSDmulti(x,
	                                   datastruc=datastruc,
									   methsMulti=methsMulti,
									   methsUni=methsUni,
									   sex.female=sex.female,
									   center=center,
									   templatevar=templatevar,
									   na.rm=na.rm,
									   ncorrection=ncorrection,
									   details=details)	  
	  out$estimates <- dimorph::logandcalcbiasMulti(out$estimates)
	}
  }
  else { #x is comparative data set, not a dimorphAds object
    uni <- NA
	{if (class(x)[1]=="data.frame" | class(x)[1]=="matrix") {
	  {if (ncol(x)==1) uni <- T
	  else if (ncol(x)>1) uni <- F}
	}
	else if (class(x)[1]=="numeric" | class(x)[1]=="integer") {
	  uni <- T
	}} # end if else
	if (is.na(uni)) stop("'x' must be a numeric vector, matrix, dataframe, or 'dimorphAds' object.")
	{if (uni) {
	  if (!is.null(struc)) {
	    struc <- struc[complete.cases(struc),,drop=F]
	    struc <- struc/struc
		if (is.na(npersample)) npersample <- nrow(struc)
	  }
	  out <- dimorph::resampleSSDuni(x,
	                                 compsex=compsex,
									 npersample=npersample,
									 nResamp=nResamp,
									 exact=exact,
									 limit=limit,
									 matchvars=matchvars,
									 replace=replace,
									 methsUni=methsUni,
									 sex.female=sex.female,
									 center=center,
									 na.rm=na.rm,
									 ncorrection=ncorrection,
									 details=details)
	  if (!is.null(struc)) out$sampleADS$struc <- struc
	  out$estimates <- dimorph::logandcalcbiasUni(out$estimates)
	}
	else if (!uni) {
	  if (is.null(datastruc)) datastruc <- "missing"
	  if (is.null(methsMulti)) methsMulti <- c("GMM")
	  out <- dimorph::resampleSSDmulti(x,
	                                   struc=struc,
									   compsex=compsex,
									   nResamp=nResamp,
									   exact=exact,
									   limit=limit,
									   matchvars=matchvars,
									   replace=replace,
									   datastruc=datastruc,
									   methsMulti=methsMulti,
									   methsUni=methsUni,
									   sex.female=sex.female,
									   center=center,
									   templatevar=templatevar,
									   na.rm=na.rm,
									   ncorrection=ncorrection,
									   details=details)	  
	  out$estimates <- dimorph::logandcalcbiasMulti(out$estimates)
	}
    } # end if else if
  }
  if (prod(is.na(out$estimates$bias))==1) out$estimates$bias <- NULL
  return(out)
}


#' @export
resampleSSDuni <- function(x,
                           compsex=NULL,
						   npersample=NA,
						   nResamp=10000,
						   exact=F,
						   limit=500000,
						   matchvars=F,
						   replace=F,
						   methsUni = c("MMR", "BDI", "BFM", "sdlog"),
						   sex.female = 1,
						   center = "geomean",
						   na.rm = T,
						   ncorrection=F,
						   details = F) {
  {if (class(x)=="dimorphAds") {
    comparative <- x$comparative
	nResamp <- x$nResamp
    struc <- x$struc
    structype <- x$strucClass
	replace <- x$replace
	exact <- x$exact
	matchvars <- x$matchvars
	compsex <- x$compsex
	sex.female <- x$sex.female
    if (!is.null(compsex)) {
      compsex <- droplevels(factor(compsex))
      x$compsex <- compsex
    }
	sampleADS <- x
	rm(x)
  }
  else {
    if (class(x)=="numeric" | class(x)=="integer") x <- data.frame(VAR=x)
    if (ncol(x) > 1) stop("This function cannot be used for multivariate datasets.")
    n <- nrow(x)
    comparative <- x
	rm(x)
	if (is.na(npersample)) npersample <- n
	if (!replace & npersample==n) warning(paste0("This combination of arguments will produce ", nResamp, "\nidentical samples equal to the full sample."))
	struc <- data.frame(x=rep(1, npersample))
	colnames(struc) <- colnames(comparative)
	sampleADS <- dimorph::getsampleaddresses(comparative=comparative,
	                                      struc=struc,
										  compsex=compsex,
										  sex.female=sex.female,
										  exact=exact,
										  limit=limit,
										  nResamp=nResamp,
										  matchvars=matchvars,
										  replace=replace)
    comparative <- sampleADS$comparative
    struc       <- sampleADS$struc
    structype   <- sampleADS$strucClass
    compsex     <- sampleADS$compsex
    sex.female  <- sampleADS$sex.female
	nResamp     <- sampleADS$nResamp
	exact       <- sampleADS$exact
    if (!is.null(compsex)) {
      compsex <- droplevels(factor(compsex))
	  sampleADS$compsex <- compsex
    }
  }}
  # check arguments
  methsUni <- match.arg(methsUni, choices = c("SSD", "MMR", "BDI", "ERM",
                                              "FMA","MoM", "BFM", 
                                              "CV", "CVsex", "sdlog", "sdlogsex"),
						several.ok=T)
  center <- match.arg(center, choices = c("geomean", "mean"))
  { # check whether compsex is present and check against comparative and methods
    if (!is.null(compsex)) {if (!length(compsex)==nrow(comparative)) stop("'compsex' must have a number of elements equal to the number of rows in 'comparative'.")}
    else if (sum(methsUni %in% c("SSD", "CVsex", "sdlogsex")) > 0) stop("Sex information must be included to use the univariate methods 'SSD', 'CVsex', or 'sdlogsex'.")
  }
  # end check of arguments
  # Declare functions that can be used with apply or lapply
  getdetails <- function(x) attr(x[[1]], "details")
  nF <- function(ads) {
    {if (is.null(compsex)) out <- NA
	else out <- sum(compsex[ads]==levels(compsex)[sex.female], na.rm=T)}
    return(out)
  } # end nF
  dimorphUniADS.data.frame <- function(ads, mUni) {
    tmp <- comparative[ads,]
    sex <- compsex
    if (!is.null(sex)) sex <- sex[ads] 
    out <- dimorph::dimorph(tmp, method = mUni, methodMulti=NA, sex=sex,
	                        #dfout=F,
	                        dfout=T,
                            sex.female=sex.female, center=center,
							na.rm=na.rm, ncorrection=ncorrection, details=F)
    return(out)
  } # end dimorphUniADS.data.frame
  # build results table
  res <- list(NULL)
  resdetails <- list(NULL)
  for (i in 1:length(methsUni)) {
    res[[i]] <- apply(sampleADS$addresses, 2, dimorphUniADS.data.frame,
		              mUni=methsUni[i])
    resdetails[[i]] <- lapply(res[[i]], FUN=getdetails)
    res[[i]] <- as.data.frame(do.call(rbind, res[[i]]))
    attr(res[[i]]$estimate,"details") <- NULL
    res[[i]] <- res[[i]][,c("estimate", "methodUni", "center", "n.specimens.overall",
		                        "proportion.female.overall", "n.specimens.realized",
								"proportion.female.realized")]
    res[[i]]$n.specimens.overall  <- nrow(sampleADS$struc)
    res[[i]]$proportion.female.overall <- apply(sampleADS$addresses, 2, nF) / res[[i]]$n.specimens.overall
    res[[i]]$subsampleID <- 1:nResamp
  }
  res <- as.data.frame(do.call(rbind, res))
  rownames(res) <- NULL
  res$methodUni <- factor(res$methodUni, levels=c("SSD", "MMR", "BDI", "ERM", "FMA", "MoM", "BFM",
                                                  "CV", "CVsex", "sdlog", "sdlogsex"))
  res$center <- factor(res$center, levels=c("geomean", "mean"))
  resdetails <- unlist(resdetails, recursive=F)
  attr(res,"details") <- resdetails
  attr(res,"estvalues") <- "raw"
  rm(resdetails)
  out <- list(estimates=res, sampleADS=sampleADS) # add in other arguments to list
  class(out) <- "dimorphResampledUni"
  attr(out, "matchvars") <- matchvars
  attr(out, "replace")   <- replace
  attr(out, "exact")     <- exact
  attr(out, "na.rm")     <- na.rm
  return(out)
}


#' @export
resampleSSDmulti <- function(x,
                             struc=NULL,
							 compsex=NULL,
							 nResamp=10000,
							 exact=F,
							 limit=500000,
							 matchvars=F,
							 replace=F,
							 datastruc="missing",
							 methsMulti = c("GMM"),
							 methsUni = c("MMR", "BDI"),
							 sex.female = 1,
							 center = "geomean",
							 templatevar = NULL,
							 na.rm = T,
							 ncorrection=F,
							 details = F) {
  {if (class(x)=="dimorphAds") {
    comparative <- x$comparative
	nResamp <- x$nResamp
    struc <- x$struc
    structype <- x$strucClass
	replace <- x$replace
	exact <- x$exact
	matchvars <- x$matchvars
	compsex <- x$compsex
	sex.female <- x$sex.female
    if (!is.null(compsex)) {
      compsex <- droplevels(factor(compsex))
      x$compsex <- compsex
    }
	sampleADS <- x
	rm(x)
  }
  else {
    comparative <- x
	rm(x)
	sampleADS <- dimorph::getsampleaddresses(comparative=comparative, 
                                            struc=struc,
											compsex=compsex,
											sex.female=sex.female,
											nResamp=nResamp,
											exact=exact,
											limit=limit,
											matchvars=matchvars,
											replace=replace)
    comparative <- sampleADS$comparative
    struc       <- sampleADS$struc
    structype   <- sampleADS$strucClass
    compsex     <- sampleADS$compsex
	sex.female  <- sampleADS$sex.female
	exact       <- sampleADS$exact
    if (!is.null(compsex)) {
      compsex <- droplevels(factor(compsex))
	  sampleADS$compsex <- compsex
    }
  }}
  # check arguments
  datastruc <- match.arg(datastruc, choices = c("missing", "complete", "both"))
  methsMulti <- match.arg(methsMulti, choices = c("GMM", "GMsize", "TM"), several.ok=T)
  #if(datastruc=="missing" & ("GMsize" %in% methsMulti)) stop("'GMsize' multivariate method cannot be performed on missing datasets.")
  if((datastruc=="missing" & ("GMsize" %in% methsMulti)) | (datastruc=="both" & ("GMsize" %in% methsMulti))) warning("'GMsize' multivariate method cannot be performed on missing datasets.")
  if(structype=="vector" & ("TM" %in% methsMulti)) warning("'TM' multivariate method cannot be performed when missing data are specified by a vector.")
  methsUni <- match.arg(methsUni, choices = c("SSD", "MMR", "BDI", "ERM",
                                              "FMA","MoM", "BFM", 
                                              "CV", "CVsex", "sdlog", "sdlogsex"),
						several.ok=T)
  center <- match.arg(center, choices = c("geomean", "mean"))
  { # check whether compsex is present and check against comparative and methods
    if (!is.null(compsex)) {if (!length(compsex)==nrow(comparative)) stop("'compsex' must have a number of elements equal to the number of rows in 'comparative'.")}
    else if (sum(methsUni %in% c("SSD", "CVsex", "sdlogsex")) > 0) stop("Sex information must be included to use the univariate methods 'SSD', 'CVsex', or 'sdlogsex'.")
  }
  # end check of arguments
  { # build table decribing complete set of methods
    if (datastruc=="both") {
      methtab <- data.frame(multi=rep(rep(methsMulti, each=length(methsUni)), 2),
                            uni=rep(rep(methsUni, length(methsMulti)), 2),
                            datstr=c(rep("complete", length(methsUni)*length(methsMulti)),
                                     rep("missing", length(methsUni)*length(methsMulti))))
    }
    else {
      methtab <- data.frame(multi=rep(methsMulti, each=length(methsUni)),
                            uni=rep(methsUni, length(methsMulti)),
                            datstr=rep(datastruc, length(methsUni)*length(methsMulti)))
    }
  }
  # remove all rows that include both GMsize and missing, since these can't happen
  methtab <- methtab[!(methtab$multi=="GMsize" & methtab$datstr=="missing"),]
  # remove all rows that include both GMM and variance-based methods, since these can't happen
  methtab <- methtab[!(methtab$multi=="GMM" & methtab$uni %in% c("CV", "CVsex", "sdlog", "sdlogsex")),]
  # if the comparative dataset has missing data, remove all rows that include complete, since these can't happen
  if (sum(is.na(comparative)) > 0) methtab <- methtab[methtab$datstr!="complete",]
  # remove all rows that include missing and one of the sex specific methods when 
  #   struc is a vector, since these can't happen
  methtab <- methtab[!(methtab$datstr=="missing" & structype=="vector" &(methtab$uni=="SSD"|methtab$uni=="CVsex"|methtab$uni=="sdlogsex")),]
  # remove all rows that include TM if structype is 'vector', since these can't happen
  if (structype=="vector") methtab <- methtab[!(methtab$multi=="TM"),]
  # reorder to put GMsize first if present
  methtab <- methtab[order(methtab$multi, decreasing=F),]
  if ("GMsize" %in% methsMulti) {
    methtabGMsize <- methtab[methtab$multi=="GMsize",]
	methtabother  <- methtab[methtab$multi!="GMsize",]
	methtab <- rbind(methtabGMsize, methtabother)
	rm(methtabGMsize, methtabother)
  }
  {if (nrow(methtab) > 0) { # At least one method combo
    rownames(methtab) <- 1:nrow(methtab)
  }
  else { # no method combos remaining
    return(NA) # Does this generate any problems for SSDtest? Doesn't appear to.
  }}
  # set up results data frame
  analyses <- paste(methtab$multi, methtab$uni, methtab$datstr, sep=".")
  #metadata <- c("n.original", "nFem.original")
  #res <- data.frame(matrix(NA, nResamp, length(metadata)+length(analyses)))
  #colnames(res) <- c(metadata, analyses)
  # Declare functions that can be used with apply or lapply
  getdetails <- function(x) attr(x[[1]], "details")
  nF <- function(ads) {
    {if (is.null(compsex)) out <- NA
	else out <- sum(compsex[ads]==levels(compsex)[sex.female], na.rm=T)}
    return(out)
  } # end nF
  dimorphMultiADS.data.frame <- function(ads, dtstr, mMulti, mUni) {
    tmp <- comparative[ads,]
    sex <- compsex
    if (!is.null(sex)) sex <- sex[ads] 
    if(dtstr=="missing") tmp <- tmp*struc
    out <- dimorph::dimorph(tmp,
	                        methodMulti = mMulti,
							method = mUni,
							sex=sex,
	                        dfout=T,
                            sex.female=sex.female,
							center=center,
							templatevar=templatevar,
							na.rm=na.rm,
							ncorrection=ncorrection,
							details=F)
    # out <- out$estimate
    return(out)
  } # end dimorphMultiADS.data.frame
  dimorphMultiADS.vector <- function(ad, dtstr, mMulti, mUni) {
    # complete method for this needs to pull together the complete set of specimens sampled 
    #  for all variables - figure out what to do if they are resampled with replacement
    if (dtstr=="complete") {
	  not.is.na <- function (x) return(!is.na(x))
      ads <- sampleADS$addresses[,ad][not.is.na(sampleADS$addresses[,ad])]
      tmp <- comparative[ads,]
      sex <- compsex
      if (!is.null(sex)) sex <- sex[ads] 
      out <- dimorph::dimorph(tmp,
                              methodMulti = mMulti,
							  method = mUni,
							  sex=sex,
							  dfout=T,
                              sex.female=sex.female,
							  center=center,
							  templatevar=templatevar,
                              na.rm=na.rm,
							  ncorrection=ncorrection,
							  details=F)
      # out <- out$estimate
    }
    if (dtstr=="missing") {
      tmp <- NULL
      for (i in names(sampleADS$adlist)) tmp[[i]] <- comparative[[i]][sampleADS$adlist[[i]][,ad]]
      nvaluesoverall <- length(unlist(tmp))
      sex <- NULL
      out <- dimorph::dimorph(tmp, 
                              methodMulti = mMulti,
							  method = mUni,
							  sex=sex,dfout=T,
							  sex.female=sex.female,
							  center=center,
							  templatevar=templatevar,
							  na.rm=na.rm,
							  ncorrection=ncorrection,
							  details = T)
	  rm(tmp)
      tmp <- NULL
      for (i in names(sampleADS$adlist)) tmp <- c(tmp, sampleADS$adlist[[i]][,ad])
      tmp <- unique(tmp)
	  nindoverall <- length(tmp)
	  nmeasurestotaloverall <- nindoverall * length(names(sampleADS$adlist))
	  nmeasuresmissingoverall <- nmeasurestotaloverall-nvaluesoverall
	  rm(tmp)
	  # calculate prop.female, proportion.missing.overall, proportion.missing.realized
	  # we can get missing overall from the specimens in adlist and ads; we can get the
	  #   prop.female and the prop missing realized by taking the variables used from
	  #   the attributes and see what specimens actually made it into the analysis
	  #varsUsed <- attr(out, "details")$vars.used
	  varsUsed <- attr(out$estimate, "details")$vars.used
      tmp <- NULL
      for (i in varsUsed) tmp <- c(tmp, sampleADS$adlist[[i]][,ad])
      tmp <- unique(tmp)
	  nind <- length(tmp)
	  nFem <- sum(compsex[tmp]==levels(compsex)[sex.female])
	  tmp2 <- NULL
      for (i in varsUsed) tmp2[[i]] <- comparative[[i]][sampleADS$adlist[[i]][,ad]]
      nvalues <- length(unlist(tmp2))
	  nmeasurestotal <- nind * length(varsUsed)
	  nmeasuresmissing <- nmeasurestotal-nvalues  
	  #attr(out, "details")$n.specimens.used <- nind
	  #attr(out, "details")$proportion.missingdata.overall <- nmeasuresmissingoverall / nmeasurestotaloverall
	  #attr(out, "details")$proportion.missingdata.realized <- nmeasuresmissing / nmeasurestotal
	  #attr(out, "details")$proportion.female <- nFem / nind
	  out$n.specimens.realized <- nind
	  out$proportion.missingdata.overall <- nmeasuresmissingoverall / nmeasurestotaloverall
	  out$proportion.missingdata.realized <- nmeasuresmissing / nmeasurestotal
	  out$proportion.female.realized <- nFem / nind
	  # out <- out$estimate
    }
    return(out)
  } # end dimorphMultiADS.vector
  { # build results table differently depending on structure type
    if (structype=="data.frame") {
	  res <- list(NULL)
	  resdetails <- list(NULL)
      for (i in 1:length(analyses)) {
        res[[i]] <- apply(sampleADS$addresses, 2, dimorphMultiADS.data.frame, dtstr=methtab$datstr[i],
                                   mMulti=methtab$multi[i], mUni=methtab$uni[i])
        resdetails[[i]] <- lapply(res[[i]], FUN=getdetails)
        res[[i]] <- as.data.frame(do.call(rbind, res[[i]]))
		attr(res[[1]]$estimate,"details") <- NULL
        res[[i]]$n.vars.overall  <- ncol(sampleADS$comparative)
        res[[i]]$n.specimens.overall  <- nrow(sampleADS$struc)
        res[[i]]$proportion.female.overall <- apply(sampleADS$addresses, 2, nF) / res[[i]]$n.specimens.overall
		res[[i]]$datastructure <- methtab$datstr[i]
		res[[i]]$subsampleID <- 1:nResamp
      }
	  names(res) <- analyses
	  res <- as.data.frame(do.call(rbind, res))
	  rownames(res) <- NULL
      res$methodUni <- factor(res$methodUni, levels=c("SSD", "MMR", "BDI", "ERM", "FMA", "MoM", "BFM",
                                                      "CV", "CVsex", "sdlog", "sdlogsex"))
      res$methodMulti <- factor(res$methodMulti, levels=c("GMsize", "GMM", "TM"))
      res$center <- factor(res$center, levels=c("geomean", "mean"))
	  res$datastructure <- factor(res$datastructure, levels=c("complete", "missing"))
	  resdetails <- unlist(resdetails, recursive=F)
      attr(res,"details") <- resdetails
	  attr(res,"estvalues") <- "raw"
	  rm(resdetails)
    } # end structype=="data.frame"
    else if (structype=="vector") {
	  res <- list(NULL)
	  resdetails <- list(NULL)
      for (i in 1:length(analyses)) {
        #res[,analyses[i]] <- apply(data.frame(x=1:nResamp), 1, dimorphMultiADS.vector, dtstr=methtab$datstr[i],
        #                           mMulti=methtab$multi[i], mUni=methtab$uni[i])
        res[[i]] <- apply(data.frame(x=1:nResamp), 1, dimorphMultiADS.vector, dtstr=methtab$datstr[i],
                                   mMulti=methtab$multi[i], mUni=methtab$uni[i])
        resdetails[[i]] <- lapply(res[[i]], FUN=getdetails)
        res[[i]] <- as.data.frame(do.call(rbind, res[[i]]))
		attr(res[[1]]$estimate,"details") <- NULL
        res[[i]]$n.vars.overall  <- ncol(sampleADS$comparative)
        #res[[i]]$n.specimens.overall  <- nrow(sampleADS$addresses) - apply(sampleADS$addresses, 2, sumNA)
        res[[i]]$n.specimens.overall  <- nrow(sampleADS$addresses) - apply(sampleADS$addresses, 2, function (x) return(sum(is.na(x))))
        res[[i]]$proportion.female.overall <- apply(sampleADS$addresses, 2, nF) / res[[i]]$n.specimens.overall
		res[[i]]$datastructure <- methtab$datstr[i]
		res[[i]]$subsampleID <- 1:nResamp
      }
	  names(res) <- analyses
	  res <- as.data.frame(do.call(rbind, res))
	  rownames(res) <- NULL
      res$methodUni <- factor(res$methodUni, levels=c("SSD", "MMR", "BDI", "ERM", "FMA", "MoM", "BFM",
                                                      "CV", "CVsex", "sdlog", "sdlogsex"))
      res$methodMulti <- factor(res$methodMulti, levels=c("GMsize", "GMM", "TM"))
      res$center <- factor(res$center, levels=c("geomean", "mean"))
	  res$datastructure <- factor(res$datastructure, levels=c("complete", "missing"))
	  resdetails <- unlist(resdetails, recursive=F)
      attr(res,"details") <- resdetails
	  attr(res,"estvalues") <- "raw"
	  rm(resdetails)
    } # end structype=="vector"
  }
  out <- list(estimates=res, sampleADS=sampleADS) # add in other arguments to list
  class(out) <- "dimorphResampledMulti"
  attr(out, "matchvars") <- matchvars
  attr(out, "replace")   <- replace
  attr(out, "exact")     <- exact
  attr(out, "datastruc") <- datastruc
  attr(out, "na.rm")     <- na.rm
  if (!("TM" %in% methtab$multi)) templatevar <- NULL
  attr(out, "templatevar") <- templatevar
  return(out)
}


#' @export
pDimorphResampledMultiOneFossil <- function(dRMcomparative, dRMfossil, fossilads=NULL, plt=T) {
  if (!(class(dRMcomparative)=="dimorphResampledMulti" & class(dRMfossil)=="dimorphResampledMulti")) {
    stop("This function requires two 'dimorphResampledMulti' objects.")
  }
  if (is.null(fossilads)) fossilads <- 1:nrow(dRMfossil$estimates)
  nComp   <- nrow(dRMcomparative$estimates)
  nFossil <- length(fossilads)
  { # adjust fossilads if the two sample sizes differ 
    if (nComp < nFossil) {
      fossilads <- sample(nFossil, nComp, replace=F) # sample without replacement fossil addresses down to comparative size
    }
    else if (nComp > nFossil) {
      nr   <- floor(nComp/nFossil)
      rd   <- nComp-nr*nFossil
      fossilads <- c(rep(fossilads, nr), sample(nFossil, rd)) # sample all an equal number of times, then sample without replacement a number of times equal to the remainder
      fossilads <- sample(fossilads, nComp) # randomize the resulting addresses
      rm(nr, rd)
    }
  } # end fossilads adjustment
  dRMfossil$estimates <- dRMfossil$estimates[fossilads,]
  estimatevars <- colnames(dRMcomparative$estimates)[colnames(dRMcomparative$estimates) %in% colnames(dRMfossil$estimates)]
  estimatevars <- estimatevars[!(estimatevars %in% c("n", "nFem"))]
  outvars <- c("p.one.sided", "p.two.sided", "mean.difference.fossil.minus.comparative", "median.difference.fossil.minus.comparative")
  d <- data.frame(matrix(NA, nComp, length(estimatevars)))
  colnames(d) <- estimatevars
  for (i in estimatevars) {
    d[,i] <- dRMfossil$estimates[,i] - dRMcomparative$estimates[,i]
  }
  out <- data.frame(matrix(NA, length(estimatevars), length(outvars)))
  colnames(out) <- outvars
  rownames(out) <- estimatevars
  pcalc <- function(x) {
    p1 <- sum(x<=0)/length(x)
    if (p1 <= 0.5) p2 <- 2*p1
    else p2 <- 2*(1-p1)
    ret <- c(p1, p2, mean(x), median(x))
    names(ret) <- c("p.one.sided", "p.two.sided",
                    "mean.difference.fossil.minus.comparative",
                    "median.difference.fossil.minus.comparative")
    return(ret)
  } # end pcalc
  out <- data.frame(t(apply(d, 2, pcalc)))
  out <- cbind(estimate=rownames(out), out)
  attr(out, "fossilads") <- fossilads
  # if 'plt'==T then plot grid of histograms
  if (plt) {
    ngridcol <- ceiling(sqrt(length(estimatevars)))
    ngridrow <- ceiling(length(estimatevars)/ngridcol)
    parOLD <- par(no.readonly=T)
    par(mfrow=c(ngridrow, ngridcol))
    nbins <- max(ceiling(nComp/100),20)
    for (i in estimatevars) {
      hist(d[,i], breaks=nbins, 
           main=paste0(i, "\nfossil - comparative taxon\np (two-sided) = ", round(out[i,"p.two.sided"],3)),
           xlab=i, cex.main=1)
      abline(v=0, col="red")
    }
    par(parOLD) 
  }
  return(out)
}

#' @export
pDimorphResampledMultiTwoFossil <- function(dRMfossil1, dRMfossil2, fossilads1=NULL, fossilads2=NULL, plt=T) {
  if (!(class(dRMfossil1)=="dimorphResampledMulti" & class(dRMfossil2)=="dimorphResampledMulti")) {
    stop("This function requires two 'dimorphResampledMulti' objects.")
  }
  if (is.null(fossilads1)) fossilads1 <- 1:nrow(dRMfossil1$estimates)
  if (is.null(fossilads2)) fossilads2 <- 1:nrow(dRMfossil2$estimates)
  nFossil1 <- length(fossilads1)
  nFossil2 <- length(fossilads2)
  { # adjust fossilads2 if the two sample sizes differ 
    if (nFossil1 < nFossil2) {
      fossilads2 <- sample(nFossil2, nFossil1, replace=F) # sample without replacement fossil addresses down to comparative size
    }
    else if (nFossil1 > nFossil2) {
      nr   <- floor(nFossil1/nFossil2)
      rd   <- nFossil1-nr*nFossil2
      fossilads2 <- c(rep(fossilads2, nr), sample(nFossil2, rd)) # sample all an equal number of times, then sample without replacement a number of times equal to the remainder
      fossilads2 <- sample(fossilads2, nFossil1) # randomize the resulting addresses
      rm(nr, rd)
    }
  } # end fossilads adjustment
  dRMfossil1$estimates <- dRMfossil1$estimates[fossilads1,]
  dRMfossil2$estimates <- dRMfossil2$estimates[fossilads2,]
  estimatevars <- colnames(dRMfossil1$estimates)[colnames(dRMfossil1$estimates) %in% colnames(dRMfossil2$estimates)]
  estimatevars <- estimatevars[!(estimatevars %in% c("n", "nFem"))]
  outvars <- c("p.one.sided", "p.two.sided", "mean.difference.fossil2.minus.fossil1", "median.difference.fossil.minus.comparative")
  d <- data.frame(matrix(NA, nFossil1, length(estimatevars)))
  colnames(d) <- estimatevars
  for (i in estimatevars) {
    d[,i] <- dRMfossil2$estimates[,i] - dRMfossil1$estimates[,i]
  }
  out <- data.frame(matrix(NA, length(estimatevars), length(outvars)))
  colnames(out) <- outvars
  rownames(out) <- estimatevars
  pcalc <- function(x) {
    p1 <- sum(x<=0)/length(x)
    if (p1 <= 0.5) p2 <- 2*p1
    else p2 <- 2*(1-p1)
    ret <- c(p1, p2, mean(x), median(x))
    names(ret) <- c("p.one.sided", "p.two.sided",
                    "mean.difference.fossil.minus.comparative",
                    "median.difference.fossil.minus.comparative")
    return(ret)
  } # end pcalc
  out <- data.frame(t(apply(d, 2, pcalc)))
  out <- cbind(estimate=rownames(out), out)
  attr(out, "fossilads1") <- fossilads1
  attr(out, "fossilads2") <- fossilads2
  # if 'plt'==T then plot grid of histograms
  if (plt) {
    ngridcol <- ceiling(sqrt(length(estimatevars)))
    ngridrow <- ceiling(length(estimatevars)/ngridcol)
    parOLD <- par(no.readonly=T)
    par(mfrow=c(ngridrow, ngridcol))
    nbins <- max(ceiling(nFossil1/100),20)
    for (i in estimatevars) {
      hist(d[,i], breaks=nbins, 
           main=paste0(i, "\nfossil2 - fossil1 taxon\np (two-sided) = ", round(out[i,"p.two.sided"],3)),
           xlab=i, cex.main=1)
      abline(v=0, col="red")
    }
    par(parOLD) 
  }
  return(out)
}

#' @export
pDimorphResampledMultiPoint <- function(dRMcomparative, dvec, plt=T) {
  if (!(class(dRMcomparative)=="dimorphResampledMulti")) {
    stop("This function requires one 'dimorphResampledMulti' object.")
  }
  nComp   <- nrow(dRMcomparative$estimates)
  estimatevars <- colnames(dRMcomparative$estimates)[colnames(dRMcomparative$estimates) %in% names(dvec)]
  estimatevars <- estimatevars[!(estimatevars %in% c("n", "nFem"))]
  outvars <- c("p.one.sided", "p.two.sided", "mean.difference.fossil.minus.comparative", "median.difference.fossil.minus.comparative")
  d <- data.frame(matrix(NA, nComp, length(estimatevars)))
  colnames(d) <- estimatevars
  for (i in estimatevars) {
    d[,i] <- dvec[i] - dRMcomparative$estimates[,i]
  }
  out <- data.frame(matrix(NA, length(estimatevars), length(outvars)))
  colnames(out) <- outvars
  rownames(out) <- estimatevars
  pcalc <- function(x) {
    p1 <- sum(x<=0)/length(x)
    if (p1 <= 0.5) p2 <- 2*p1
    else p2 <- 2*(1-p1)
    ret <- c(p1, p2, mean(x), median(x))
    names(ret) <- c("p.one.sided", "p.two.sided",
                    "mean.difference.fossil.minus.comparative",
                    "median.difference.fossil.minus.comparative")
    return(ret)
  } # end pcalc
  out <- data.frame(t(apply(d, 2, pcalc)))
  out <- cbind(estimate=rownames(out), out)
  # if 'plt'==T then plot grid of histograms
  if (plt) {
    ngridcol <- ceiling(sqrt(length(estimatevars)))
    ngridrow <- ceiling(length(estimatevars)/ngridcol)
    parOLD <- par(no.readonly=T)
    par(mfrow=c(ngridrow, ngridcol))
    nbins <- max(ceiling(nComp/100),20)
    for (i in estimatevars) {
      hist(d[,i], breaks=nbins, 
           main=paste0(i, "\nfossil - comparative taxon\np (two-sided) = ", round(out[i,"p.two.sided"],3)),
           xlab=i, cex.main=1)
      abline(v=0, col="red")
    }
    par(parOLD) 
  }
  return(out)
}


#' @export
dimorphADS <- function(ads, x, method="SSD", methodMulti="GMM", sex=NULL, sex.female=1,
                       center="geomean", templatevar=NULL, na.rm=T, details=F, dfout=F) {
  return(dimorph(x=x, method=method, methodMulti=methodMulti, sex=sex, sex.female=sex.female,
                 center=center, templatevar=templatevar, ads=ads, na.rm=na.rm, details=details, dfout=dfout))
}


#' @export
build_dimorphAds <- function(addresses, comparative, compsex, struc=NULL, replace) {
  nResamp <- ncol(addresses)
  adlist <- NULL
  if (is.null(struc)) struc <- data.frame(matrix(1, nrow(addresses), ncol(comparative)))
  colnames(struc) <- colnames(comparative)
  strucClass <- "data.frame"
  out <- list(nResamp=nResamp,
              addresses=addresses,
              adlist=adlist,
              comparative=comparative,
              compsex=compsex,
              replace=replace,
              struc=struc,
              strucClass=strucClass)
  class(out) <- "dimorphAds"
  return(out)
}

#' @export
build_dimorphAdsStruc <- function(struc=NULL, addresses, comparative, compsex, replace) {
  nResamp <- ncol(addresses)
  adlist <- NULL
  if (is.null(struc)) struc <- data.frame(matrix(1, nrow(addresses), ncol(comparative)))
  colnames(struc) <- colnames(comparative)
  strucClass <- "data.frame"
  out <- list(nResamp=nResamp,
              addresses=addresses,
              adlist=adlist,
              comparative=comparative,
              compsex=compsex,
              replace=replace,
              struc=struc,
              strucClass=strucClass)
  class(out) <- "dimorphAds"
  return(out)
}

#' @export
logandcalcbiasUni <- function(resDF) {
  resDF$bias <- NA
  if (attr(resDF, "estvalues")=="raw") {
    resDF$estimate[resDF$methodUni %in% c("SSD", "MMR", "BDI", "MoM", "FMA", "BFM", "ERM")] <- 
       log(resDF$estimate[resDF$methodUni %in% c("SSD", "MMR", "BDI", "MoM", "FMA", "BFM", "ERM")])
    attr(resDF, "estvalues") <- "logged"
  }
  for (i in levels(resDF$methodUni)) {
    for (j in levels(resDF$center)) {
      if (i %in% c("MMR", "BDI", "MoM", "FMA", "BFM", "ERM")) {
	    actual <- resDF[resDF$methodUni=="SSD" & resDF$center==j, "estimate"]
		if (length(actual > 0)) resDF[resDF$methodUni==i & resDF$center==j, "bias"] <- resDF[resDF$methodUni==i & resDF$center==j, "estimate"] - actual
        rm(actual)
	  }
	  else if (i=="CV") {
	    actual <- resDF[resDF$methodUni=="CVsex" & resDF$center=="mean", "estimate"]
	    if (length(actual > 0)) resDF[resDF$methodUni==i & resDF$center==j, "bias"] <- resDF[resDF$methodUni==i & resDF$center==j, "estimate"] - actual
	  }
	  else if (i=="sdlog") {
	    actual <- resDF[resDF$methodUni=="sdlogsex" & resDF$center=="geomean", "estimate"]
	    if (length(actual > 0)) resDF[resDF$methodUni==i & resDF$center==j, "bias"] <- resDF[resDF$methodUni==i & resDF$center==j, "estimate"] - actual
	  }
	}
  }
  return(resDF)
}

#' @export
logandcalcbiasMulti <- function(resDF) {
  resDF$bias <- NA
  if (class(resDF$methodUni)=="character") resDF$methodUni <- factor(resDF$methodUni,
                levels=c("SSD", "MMR", "BDI", "ERM", "MoM", "FMA", "BFM", "CV", "CVsex", "sdlog", "sdlogsex"))
  if (class(resDF$center)=="character") resDF$center <- factor(resDF$center, levels=c("geomean", "mean"))
  if (class(resDF$methodMulti)=="character") resDF$methodMulti <- factor(resDF$methodMulti,
                                                                         levels=c("GMsize", "GMM", "TM"))
  if (class(resDF$datastructure)=="character") resDF$datastructure <- factor(resDF$datastructure,
                                                                             levels=c("complete", "missing"))
  if (attr(resDF, "estvalues")=="raw") {
    resDF$estimate[resDF$methodUni %in% c("SSD", "MMR", "BDI", "ERM", "MoM", "FMA", "BFM")] <- 
       log(resDF$estimate[resDF$methodUni %in% c("SSD", "MMR", "BDI", "ERM", "MoM", "FMA", "BFM")])
    attr(resDF, "estvalues") <- "logged"
  }
  for (y in levels(resDF$datastructure)) {
    for (z in levels(resDF$methodMulti)) {
      for (i in levels(resDF$methodUni)) {
        for (j in levels(resDF$center)) {
          if (i %in% c("MMR", "BDI", "ERM", "MoM", "FMA", "BFM")) {
    	    actual <- resDF[resDF$methodUni=="SSD" & resDF$center==j & 
			                resDF$datastructure==y & resDF$methodMulti==z, "estimate"]
			if (length(actual > 0)) resDF[resDF$methodUni==i & resDF$center==j & 
			                        resDF$datastructure==y & resDF$methodMulti==z, "bias"] <- 
			                          resDF[resDF$methodUni==i & resDF$center==j & 
			                          resDF$datastructure==y & resDF$methodMulti==z, "estimate"] - actual
            rm(actual)
	      }
	      else if (i=="CV") {
	        actual <- resDF[resDF$methodUni=="CVsex" & resDF$center=="mean" & 
			                resDF$datastructure==y & resDF$methodMulti==z, "estimate"]
	        if (length(actual > 0)) resDF[resDF$methodUni==i & resDF$center==j & 
			                        resDF$datastructure==y & resDF$methodMulti==z, "bias"] <- 
			                          resDF[resDF$methodUni==i & resDF$center==j & 
			                          resDF$datastructure==y & resDF$methodMulti==z, "estimate"] - actual
	      }
	      else if (i=="sdlog") {
	        actual <- resDF[resDF$methodUni=="sdlogsex" & resDF$center=="geomean" & 
			                resDF$datastructure==y & resDF$methodMulti==z, "estimate"]
	        if (length(actual > 0)) resDF[resDF$methodUni==i & resDF$center==j & 
			                        resDF$datastructure==y & resDF$methodMulti==z, "bias"] <- 
			                          resDF[resDF$methodUni==i & resDF$center==j & 
			                          resDF$datastructure==y & resDF$methodMulti==z, "estimate"] - actual
	      }
	    } # end j
      } # end i
    } # end z
  } # end y
  return(resDF)
}

#' @export
# Get confidence intervals for full ape data set
getCI <- function(x, conf.level=0.95, alternative="two.sided", na.rm=T) {
  if (!na.rm & sum(is.na(x)) > 0) {
    warning("Resampled values contain NAs.  Set 'na.rm' to TRUE to calculate confidence interval.")
	out <- c(lower_lim=NA, upper_lim=NA, proportionNA=NA)
	out <- as.numeric(out)
	return(out)
  }
  nOrig <- length(x)
  x <- sort(x)
  alternative <- match.arg(alternative, choices=c("two.sided", "less", "greater"))
  alpha=1-conf.level
  n=length(x)
  if (n==0) {
	out <- c(lower_lim=NA, upper_lim=NA, proportionNA=1)
	out <- as.numeric(out)
	return(out) 
  }
  proportionNA <- 1-(n/nOrig)
  ll=1
  ul=n
  if (alternative=="two.sided") {
    chuck <- floor(alpha/2*n)
	ll <- chuck+1
	ul <- n-chuck
	out <- c(lower_lim=x[ll], upper_lim=x[ul], proportionNA=proportionNA)
  }
  else if (alternative=="less") {
    chuck <- floor(alpha*n)
	ul <- n-chuck
	out <- c(lower_lim=-Inf, upper_lim=x[ul], proportionNA=proportionNA)
  }
  else {
    chuck <- floor(alpha*n)
	ll <- chuck+1
	out <- c(lower_lim=x[ll], upper_lim=Inf, proportionNA=proportionNA)
 }
 return(out)
}

