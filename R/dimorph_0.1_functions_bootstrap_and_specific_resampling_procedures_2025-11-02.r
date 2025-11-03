#' Bootstrap dimorphism confidence intervals in univariate or multivariate sample
#' 
#' Function to generate confidence intervals for actual or estimated dimorphism in a univariate or multivariate sample
#'   by percentile bootstrapping; i.e., sampling with replacement multiple times to a sample size equal to the original 
#'   dataset sample size, then discarding the highest and lowest (alpha/2) x 100 percent of resampled values to 
#'   find confidence limits (or discard the highest or the lowest alpha x 100 percent of values in the case of 
#'   one-sided confidence intervals).
#' @param x A matrix or data frame of measurements from a comparative sample, with rows
#'   corresponding to individual specimens and columns corresponding to size variables.  Sex data
#'   should not be included.
#' @param sex A vector indicating sex for the individuals in \code{comparative}.  Defaults to \code{NULL}.
#' @param methsUni A character vector specifying the univariate method(s) used to calculate or estimate dimorphism. 
#'   See \code{\link{dimorph}} for options.
#' @param methsMulti A character vector specifying the multivariate method(s) used to 
#'   calculate or estimate dimorphism.  Note that regardless of the value of this argument,
#'   multivariate estimation procedures will only be carried out if \code{x} is a multivariate 
#'   dataset.  See \code{\link{dimorph}} for options.
#' @param nResamp Integer specifying the number of resampling iterations if Monte Carlo sampling is used.
#' @param exact Logical scalar specifying whether to sample all unique combinations of resampled datasets at the same
#'   sample size as \code{x} from \code{x} when sampling with replacement.  Defaults to \code{FALSE}.  If set to 
#'   \code{FALSE}, or if set to \code{TRUE} and the number of unique combinations exceeds \code{limit}, then Monte  
#'   Carlo sampling is used instead.
#' @param limit Integer setting the upper limit on the number of unique combinations allowable 
#'   for exact resampling.  If exact resampling would produce more resampled datasets than this number, 
#'   Monte Carlo resampling is used instead.  Defaults to 50,000.
#' @param center A character string specifying the method used to calculate a mean, either \code{"geomean"} 
#'   (default) which uses the geometric mean, or \code{"mean"} which uses the arithmetic mean.  More broadly, 
#'   \code{"geomean"} indicates analyses are conducted in logarithmic data space and \code{"mean"} indicates 
#'   analyses are conducted in raw data space.  Some methods can only be applied in one domain or the other: 
#'   \code{"CV"} and \code{"CVsex"} are always calculated in raw data space and \code{center} will be set to 
#'   \code{"mean"} for these methods regardless of the value set by the user; \code{"MoM"}, \code{"sdlog"}, 
#'   and \code{"sdlogsex"} are always calculated in logarithmic data space and \code{center} will be set to 
#'   \code{"geomean"} for these methods regardless of the value set by the user.
#' @param sex.female An integer scalar (1 or 2) specifying which level of \code{sex} 
#'   corresponds to female.  Ignored if \code{sex} is \code{NULL}.  Defaults to 1.
#' @param na.rm A logical scalar indicating whether \code{NA} values should be stripped before
#'   the computation proceeds in univariate analyses.  Not relevant for multivariate analyses.  Defaults to \code{TRUE}.
#' @param ncorrection A logical scalar indicating whether to apply Sokal and Braumann's (1980) 
#'   size correction factor to CV estimates.  Defaults to \code{FALSE}.
#' @param struc Not generally relevant for users, this argument is sometimes applicable when \code{bootdimorph}
#'   is called by other functions.  Defaults to \code{NULL}.  See \code{\link{resampleSSD}} for more details. 
#' @param datastruc Not generally relevant for users, this argument is sometimes applicable when \code{bootdimorph}
#'   is called by other functions.  Defaults to \code{NULL}.  See \code{\link{resampleSSD}} for more details. 
#' @param templatevar A character object or integer value specifying the name or column number 
#'   of the variable in \code{x} to be estimated using the template method.  Ignored if template method 
#'   is not used.  Defaults to \code{NULL}.
#' @param alternative A character object specifying whether to calculate two-sided (\code{"two.sided"}), or 
#'   one-sided (\code{"less"} or \code{"greater"}) confidence intervals.  Defaults to \code{"two.sided"}.
#' @param conf.level Value between zero and one setting the confidence level, equal to 1-alpha. Defaults
#'   to \code{0.95}.
#' @return A list of class \code{dimorphResampledUni} or \code{dimorphResampledMulti} containing a dataframe
#'   with resampled dimorphism estimates and a \code{dimorphAds} object containg resampled addresses produced
#'   by \code{\link{getsampleaddresses}}.  Printing this object provides confidence intervals for all estimators calculated,
#'   and confidence intervals for bias from sample SSD if relevant.  Plotting this object produces violin plots for
#'   all bootstrapped distributions and lines indicating confidence limits.
#' @seealso \code{\link{centers}}, \code{\link{confint}}, \code{\link{resampleSSD}}, 
#'   \code{\link{dimorph}}, \code{\link{getsampleaddresses}}
#' @examples
#' gor <- apelimbart[apelimbart$Species=="Gorilla gorilla",]
#' hom <- apelimbart[apelimbart$Species=="Homo sapiens",]
#' pan <- apelimbart[apelimbart$Species=="Pan troglodytes",]
#' hyl <- apelimbart[apelimbart$Species=="Hylobates lar",]
#' outcomeUgor <- bootdimorph(gor[,"FHSI", drop=FALSE], sex=gor$Sex, nResamp=100)
#' outcomeUhom <- bootdimorph(hom[,"FHSI", drop=FALSE], sex=hom$Sex, nResamp=100)
#' outcomeUpan <- bootdimorph(pan[,"FHSI", drop=FALSE], sex=pan$Sex, nResamp=100)
#' outcomeUhyl <- bootdimorph(hyl[,"FHSI", drop=FALSE], sex=hyl$Sex, nResamp=100)
#' outcomeUgor
#' confint(outcomeUgor)
#' confint(outcomeUgor, conf.level=0.8, alternative="greater")
#' plot(outcomeUgor)
#' plot(outcomeUgor, exclude="FMA") # exclude one or more methods from plot
#' outcomeUhom
#' plot(outcomeUhom)
#' outcomeUpan
#' plot(outcomeUpan)
#' outcomeUhyl
#' plot(outcomeUhyl)
#' confint(outcomeUgor, type="bias")
#' plot(outcomeUgor, type="bias")
#' plot(outcomeUhom, type="bias")
#' plot(outcomeUpan, type="bias")
#' plot(outcomeUhyl, type="bias")
#' @export
bootdimorph <- function(x, 
						sex=NULL,
						methsUni=c("SSD", "MMR", "BDI", "ERM", "FMA", "MoM", "BFM",
                                   "CV", "CVsex", "sdlog", "sdlogsex"),
                        methsMulti="GMM",
                        nResamp=1000,
						exact=F,
						limit=50000,
                        center="geomean",
						sex.female=1,
						na.rm=T,
						ncorrection=F,
						struc=NULL,
						datastruc=NULL,
						templatevar=NULL,
						alternative="two.sided",
						conf.level=0.95) { 
  alternative <- match.arg(alternative, choices=c("two.sided", "less", "greater"))
  methsUni <- match.arg(methsUni, choices=c("SSD", "MMR", "BDI", "ERM", "FMA", "MoM", "BFM",
                                            "CV", "CVsex", "sdlog", "sdlogsex"), several.ok=T)
  methsMulti <- match.arg(methsMulti, choices=c("GMsize", "GMM", "TM"), several.ok=T)
  center <- match.arg(center, choices=c("geomean", "mean"))
  if (class(conf.level)!="numeric") stop("'conf.level' must be a number between 0 and 1.")
  if (conf.level <= 0 | conf.level >= 1) stop("'conf.level' must be a number between 0 and 1.")
  # code to check that x is either a vector of values or a data.frame or matrix of values
  uni <- NA
  {if (class(x)[1]=="data.frame" | class(x)[1]=="matrix") {
    {if (ncol(x)==1) {
	  uni <- T
	  if (!("data.frame" %in% class(x))) {
	    if (is.null(colnames(x))) colnames(x) <- "VAR"
	    x <- data.frame(x)
	  }
	  # check for missing data, proceed according to na.rm, generate warning or error
	  if (sum(is.na(x)) > 0 | sum(is.na(sex)) > 0) {
	    if (!na.rm) stop("There are missing data and 'na.rm' is FALSE.")
	    else {
	      rmads <- unique(c(which(is.na(x[,1])), which(is.na(sex))))
		  x <- x[-rmads, , drop=FALSE]
		  if (!is.null(sex)) {
		    sex <- sex[-rmads]
		    sex <- droplevels(factor(sex))
		  }
		  warning(paste0(length(rmads), " observations have been removed due to missing data."))
	    }
	  }
	  #xnames <- rownames(x)
	  #x <- unlist(x)
	  #names(x) <- xnames
	  #rm(xnames)
	}
    else if (ncol(x)>1) uni <- F}
  }
  else if (class(x)[1]=="numeric" | class(x)[1]=="integer") {
    x <- data.frame(VAR=x)
    uni <- T
	# check for missing data, proceed according to na.rm, generate warning or error
	if (sum(is.na(x)) > 0 | sum(is.na(sex)) > 0) {
	  if (!na.rm) stop("There are missing data and 'na.rm' is FALSE.")
	  else {
	    rmads <- unique(c(which(is.na(x[,1])), which(is.na(sex))))
		x <- x[-rmads, , drop=FALSE]
		if (!is.null(sex)) {
		  sex <- sex[-rmads]
		  sex <- droplevels(factor(sex))
		}
		warning(paste0(length(rmads), " observations have been removed due to missing data."))
	  }
	}
  }
  else if (class(x)[1]=="list") uni <- F}
  if (is.na(uni)) stop("'x' must be a numeric vector, matrix, dataframe, or list.")
  if (is.null(sex)) {
	methsUni <- methsUni[!(methsUni %in% c("SSD", "CVsex", "sdlogsex"))]
	if (length(methsUni)==0) stop("Methods 'SSD', 'CVsex', and 'sdlogsex' require 'sex' data.")
  }
  # get resampled values
  {if (uni) {
    methsMulti <- NULL
	if (na.rm) x <- x[stats::complete.cases(x),,drop=F]
	res <- dimorph:::resampleSSDuni(x=x,
	                      compsex=sex,
						  npersample=NA,
	                      nResamp=nResamp,
						  exact=exact,
						  limit=limit,
						  matchvars=F,
						  replace=T,
						  methsUni=methsUni,
						  sex.female=sex.female,
						  center=center,
						  na.rm=na.rm,
						  ncorrection=ncorrection,
						  details=F)
    attributes(res$estimates)$details <- NULL
    attr(res$estimates, "estvalues") <- "raw"
    res$estimates <- dimorph:::logandcalcbiasUni(res$estimates)
    CI <- tapply(res$estimates$estimate, INDEX=res$estimates$methodUni:res$estimates$center,
	             FUN=dimorph:::getCI, alternative=alternative, conf.level=conf.level)
	CI <- do.call(rbind, CI)
	CImeta <- do.call(rbind, strsplit(rownames(CI), ":"))
	CI <- data.frame(methodUni=CImeta[,1], center=CImeta[,2],
	                 lower_lim=CI[,1], upper_lim=CI[,2])
	CI$methodUni <- factor(CI$methodUni, levels=levels(res$estimates$methodUni))
	CI$center <- factor(CI$center, levels=levels(res$estimates$center))
	rownames(CI) <- NULL
	res$CI <- CI
	if (sum(is.na(res$estimates$bias)) != length(res$estimates$bias)) {
      CIbias <- tapply(res$estimates$bias, INDEX=res$estimates$methodUni:res$estimates$center,
	                   FUN=dimorph:::getCI, alternative=alternative, conf.level=conf.level)
	  CIbias <- do.call(rbind, CIbias)
	  CImeta <- do.call(rbind, strsplit(rownames(CIbias), ":"))
	  CIbias <- data.frame(methodUni=CImeta[,1], center=CImeta[,2],
	                 lower_lim=CIbias[,1], upper_lim=CIbias[,2])
	  CIbias$methodUni <- factor(CIbias$methodUni, levels=levels(res$estimates$methodUni))
	  CIbias$center <- factor(CIbias$center, levels=levels(res$estimates$center))
	  CIbias <- CIbias[!(CIbias$methodUni %in% c("SSD", "CVsex", "sdlogsex")),]
	  rownames(CIbias) <- NULL
	  res$CIbias <- CIbias
	}
	attr(res, "resampling") <- "bootstrap"
	attr(res, "conf.level") <- conf.level
	attr(res, "alternative") <- alternative
    return(res)
  }
  else if (!uni) {
    if (is.null(struc)) {
	  struc <- x/x
	  datastruc <- "complete"
	}
    datastruc <- match.arg(datastruc, choices=c("complete", "missing", "both"))
	res <- dimorph:::resampleSSDmulti(x=x,
	                        struc=struc,
							compsex=sex,
	                        nResamp=nResamp,
							exact=exact,
							limit=limit,
							matchvars=F,
							replace=T,
							methsUni=methsUni,
							methsMulti=methsMulti,
							datastruc=datastruc,
							templatevar=templatevar,
							sex.female=sex.female,
							center=center,
							na.rm=na.rm,
							ncorrection=ncorrection,
							details=F)
    attributes(res$estimates)$details <- NULL
    attr(res$estimates, "estvalues") <- "raw"
    res$estimates <- dimorph:::logandcalcbiasMulti(res$estimates)
    CI <- tapply(res$estimates$estimate,
	             INDEX=res$estimates$methodUni:res$estimates$methodMulti:res$estimates$center:res$estimates$datastructure,
				 FUN=dimorph:::getCI, alternative=alternative, conf.level=conf.level)
	CI <- do.call(rbind, CI)
	CImeta <- do.call(rbind, strsplit(rownames(CI), ":"))
	CI <- data.frame(methodUni=CImeta[,1], methodMulti=CImeta[,2], center=CImeta[,3],
	                 datastructure=CImeta[,4], proportionNA=CI[,3],
					 lower_lim=CI[,1], upper_lim=CI[,2])
	CI$methodUni <- factor(CI$methodUni, levels=levels(res$estimates$methodUni))
	CI$methodMulti <- factor(CI$methodMulti, levels=levels(res$estimates$methodMulti))
	CI$center <- factor(CI$center, levels=levels(res$estimates$center))
	CI$datastructure <- factor(CI$datastructure, levels=levels(res$estimates$datastructure))
	rownames(CI) <- NULL
	res$CI <- CI
	if (sum(is.na(res$estimates$bias)) != length(res$estimates$bias)) {
      CIbias <- tapply(res$estimates$bias,
	                   INDEX=res$estimates$methodUni:res$estimates$methodMulti:res$estimates$center:res$estimates$datastructure,
					   FUN=dimorph:::getCI, alternative=alternative, conf.level=conf.level)
	  CIbias <- do.call(rbind, CIbias)
	  CImeta <- do.call(rbind, strsplit(rownames(CIbias), ":"))
	  CIbias <- data.frame(methodUni=CImeta[,1], methodMulti=CImeta[,2], center=CImeta[,3],
	                       datastructure=CImeta[,4], proportionNA=CIbias[,3],
						   lower_lim=CIbias[,1], upper_lim=CIbias[,2])
	  CIbias$methodUni <- factor(CIbias$methodUni, levels=levels(res$estimates$methodUni))
	  CIbias$methodMulti <- factor(CIbias$methodMulti, levels=levels(res$estimates$methodMulti))
	  CIbias$center <- factor(CIbias$center, levels=levels(res$estimates$center))
	  CIbias$datastructure <- factor(CIbias$datastructure, levels=levels(res$estimates$datastructure))
	  CIbias <- CIbias[!(CIbias$methodUni %in% c("SSD", "CVsex", "sdlogsex")),]
	  rownames(CIbias) <- NULL
	  res$CIbias <- CIbias
	}
	attr(res, "resampling") <- "bootstrap"
	attr(res, "conf.level") <- conf.level
	attr(res, "alternative") <- alternative
    return(res)    
  }}
}

#' @noRd
bootsampleaddresses <- function(x,
                                nResamp=NA) {
  # don't allow exact resampling for now - the number of resampled values for the smallest fossil dataset
  #   (A. afarensis downsampled to minimum of afarensis and africanus sample sizes) is over 9 million.
  #   Eventually code this in, but don't worry now, especially since I'll need to adjust the exact list and 
  #   count if there was sampling with replacement in the first case.
  exact <- F
  # move 'exact' to arguments list later
  if (!inherits(x, "dimorphAds")) stop("'x' must be a 'dimorphAds' object.")
  #if (class(x)!="dimorphAds") stop("'x' must be a 'dimorphAds' object.")
  x$originalstrucClass <- x$strucClass 
  nvar <- ncol(x$comparative)
  if (x$strucClass=="data.frame") { # generates adlist
    #nspeccomp <- nrow(x$comparative)
    #nspecstruc <- nrow(x$struc)
	#nvar <- ncol(x$comparative)
	nspecpervar <- apply(x$struc, 2, function(x) sum(!is.na(x)))
	# if original set of sample addresses was done without replacement, then
	#nexactpersample <- choose(2*nspecpervar-1, nspecpervar) # combination WITH replacement (n+k-1, choose k)
	x$adlist <- list(NULL)
	for (i in 1:nvar) {
	  goodAds <- which(!is.na(x$struc[,i]))
	  x$adlist[[i]] <- x$addresses[goodAds,,drop=F]
	}
	names(x$adlist) <- colnames(x$comparative)
	strucrownames <- rownames(x$struc)
	x$struc <- nspecpervar
	x$strucClass <- "vector"
  } # end data.frame procedure
  # from here on data.frame and vector procedure are the same
  if (!is.na(nResamp)) {
    if(nResamp > x$nResamp) {
    npersubsample <- ceiling(nResamp/x$nResamp)
	for (i in 1:nvar) {
	  x$adlist[[i]] <- matrix(rep(x$adlist[[i]], npersubsample), nrow(x$adlist[[i]]),
	                              ncol(x$adlist[[i]])*npersubsample, byrow=F)
	}    
  }}
  for (i in 1:nvar) { # resample with replacement in each subsample
    x$adlist[[i]] <- apply(x$adlist[[i]], 2, sample, size=nrow(x$adlist[[i]]), replace=T)
  }
  x$nResamp <- ncol(x$adlist[[1]])
  if (!is.na(nResamp)) {
  if (x$nResamp > nResamp) {
    keepads <- sample(1:x$nResamp, size=nResamp, replace=F)
	for (i in 1:nvar) { # resample with replacement in each subsample
	  x$adlist[[i]] <- x$adlist[[i]][,keepads]
	}
	x$nResamp <- nResamp
  }} 
  { # start if else
  if (nvar==1) { # if univariate
    x$addresses <- x$adlist[[1]]
	x$adlist <- NULL
	x$struc <- data.frame(x=rep(1, x$struc))
	colnames(x$struc) <- colnames(x$comparative)
	if (x$originalstrucClass=="data.frame") rownames(x$struc) <- strucrownames
	x$strucClass <- "data.frame"
  } # end univariate
  
  else { # if multivariate
    nmax <- sum(x$struc)
	#getlistreplaceF <- function(i) {
	#  tmp <- NULL
	#  for (j in 1:length(x$adlist)) tmp <- c(tmp, x$adlist[[j]][,i])
	#  tmp <- unique(tmp)
	#  ret <- rep(NA, nmax)
	#  ret[1:length(tmp)] <- tmp
	#  return(ret)
	#}
	getlistreplaceT <- function(i) {
	  tmp <- NULL
	  for (j in 1:length(x$adlist)) {
	    tmp <- c(tmp, table(x$adlist[[j]][,i]))
	  }
	  tmp <- sort(tmp, decreasing=T)
	  tmp <- tmp[unique(names(tmp))]
	  tmp2 <- NULL
	  for (k in 1:length(tmp)) tmp2 <- c(tmp2, rep(as.integer(names(tmp)[k]), tmp[k]))
	  ret <- rep(NA, nmax)
	  ret[1:length(tmp2)] <- tmp2
	  return(ret)
	}
	#if (!replace) ads <- apply(data.frame(x=1:x$nResamp), 1, getlistreplaceF)
	#if (replace)  ads <- apply(data.frame(x=1:x$nResamp), 1, getlistreplaceT)
	x$addresses <- apply(data.frame(x=1:x$nResamp), 1, getlistreplaceT)
	x$addresses <- x$addresses[apply(x$addresses, 1, function (x) return(sum(is.na(x)))) < x$nResamp,] # remove completely empty rows
  } # end multivariate
  } # end if else
  # add more tags to dimorphAds object to indicate that it's been resampled again
  x$rebootstrapped <- T
  return(x)
}

#' Significance Tests for Sexual Size Dimorphism Estimates
#' 
#' Function to calculate \emph{p}-values for pairwise comparisons of dimorphism estimates for two or more univariate 
#'   or multivariate samples.  The first step in this process identifies the minimum sample size (and missing data 
#'   structure) present in the samples of \code{fossil}; for multivariate datasets with missing data, datasets are 
#'   restricted to individuals with measurements in variables present in all samples.  All datasets are then resampled to
#'   the sample size (and missing data structure) of the minimum sample.  This function generates both
#'   two-sided and one-sided \emph{p}-values for each pair of samples.
#' @param fossil A list of matrices or data frames of measurements from fossil sample(s) (or other test samples) to be
#'   compared in a series of pairwise resampling tests, with rows corresponding to individual specimens and 
#'   columns corresponding to size variables.  Sex data should not be included.  Some data can be missing. 
#'   When a single fossil sample is present, by default null distributions of resampled values from comparative data sets 
#'   will be compared to a single point estimate generated for the fossil sample (one point estimate per estimation 
#'   method).  \code{fossil} may be \code{NULL}, in which case all samples should be provided in \code{comp}, and all 
#'   samples will be resampled.
#' @param comp A list of matrices or data frames of measurements from comparative sample(s) to be
#'   compared in a series of pairwise resampling tests, with rows corresponding to individual specimens and 
#'   columns corresponding to size variables.  Sex data should not be included.  All data sets must be complete
#'   for all measurements.
#' @param fossilsex A list of vectors indicating sex for the individuals in each of the samples in \code{fossil}. 
#'   Defaults to \code{NULL}.
#' @param compsex A list of vectors indicating sex for the individuals in each of the samples in \code{comp}. 
#'   Defaults to \code{NULL}.
#' @param methsUni A character vector specifying the univariate method(s) used to calculate or estimate dimorphism. 
#'   See \code{\link{dimorph}} for options.  Defaults to \code{c("MMR", "BDI")}.
#' @param methsMulti A character vector specifying the multivariate method(s) used to 
#'   calculate or estimate dimorphism.  Note that regardless of the value of this argument,
#'   multivariate estimation procedures will only be carried out if \code{fossil} and \code{comp} are multivariate 
#'   datasets.  See \code{\link{dimorph}} for options.  Defaults to \code{"GMM"}.
#' @param replace Logical scalar specifying whether to sample from comparative datasets with replacement or not.
#'   Defaults to \code{FALSE}. 
#' @param rebootstrap Logical scalar specifying whether to add an additional step after initial resampling in 
#'   which resampled addresses for both comparative and fossil datasets are bootstrapped (i.e., sampled with 
#'   replacement to an equal sample sizes as the initial set of resampled addresses).  This procedure implements 
#'   the highly conservative "resampled extinct distribution method" of Gordon et al. (2008). Defaults to 
#'   \code{FALSE}.  Warning: setting \code{rebootstrap} to \code{TRUE} will drastically reduce power!
#' @param fullsamplesboot Logical scalar.  If all samples are complete (no missing data) and \code{fullsamplesboot} is 
#'   set to \code{TRUE}, rather than downsampling to the minimum sample in \code{fossil}, all samples in \code{fossil} are
#'   and \code{comp} are bootstrapped (i.e., sampled with replacement at their full sample size).  Defaults to 
#'   \code{FALSE}.  Note that for univariate analyses, \code{NA}s are removed by default so all datasets (including 
#'   fossils) are considered complete.  Therefore setting \code{fullsamplesboot} to \code{TRUE} for univariate analyses 
#'   with fossils will bootstrap the fossil dataset rather than generating a single point estimate from the fossil dataset.
#' @param nResamp Integer specifying the number of resampling iterations 
#'   if Monte Carlo sampling is used.  Defaults to \code{1,000}.
#' @param exactcomp Logical scalar specifying for samples in \code{comp} whether to sample all possible unique 
#'   combinations of resampled datasets for the minimum data structure present in \code{fossil}. 
#'   If set to \code{FALSE}, or if set to \code{TRUE} and the number of unique combinations exceeds \code{limit}, 
#'   then Monte Carlo sampling is used instead.  Defaults to \code{TRUE}.
#' @param exactfossil Logical scalar specifying for samples in \code{fossil} whether to sample all 
#'   possible unique combinations of resampled datasets for the minimum data structure present in \code{fossil}. 
#'   If set to \code{FALSE}, or if set to \code{TRUE} and the number of unique combinations exceeds \code{limit}, 
#'   then Monte Carlo sampling is used instead.  Defaults to \code{TRUE}.
#' @param limit Integer setting the upper limit on the number of unique combinations allowable 
#'   for exact resampling.  If exact resampling would produce more resampled datasets than this number, 
#'   Monte Carlo resampling is used instead.  Defaults to \code{50,000}.
#' @param center A character string specifying the method used to calculate a mean, either \code{"geomean"} 
#'   (default) which uses the geometric mean, or \code{"mean"} which uses the arithmetic mean.  More broadly, 
#'   \code{"geomean"} indicates analyses are conducted in logarithmic data space and \code{"mean"} indicates 
#'   analyses are conducted in raw data space.  Some methods can only be applied in one domain or the other: 
#'   \code{"CV"} and \code{"CVsex"} are always calculated in raw data space and \code{center} will be set to 
#'   \code{"mean"} for these methods regardless of the value set by the user; \code{"MoM"}, \code{"sdlog"}, 
#'   and \code{"sdlogsex"} are always calculated in logarithmic data space and \code{center} will be set to 
#'   \code{"geomean"} for these methods regardless of the value set by the user.
#' @param sex.female An integer scalar (1 or 2) specifying which level of \code{sex} 
#'   corresponds to female.  Ignored if \code{sex} is \code{NULL}.  Defaults to 1.
#' @param na.rm A logical scalar indicating whether \code{NA} values should be stripped before
#'   the computation proceeds in univariate analyses.  Not relevant for multivariate analyses.  Defaults to \code{TRUE}.
#' @param ncorrection A logical scalar indicating whether to apply Sokal and Braumann's (1980) 
#'   size correction factor to CV estimates.  Defaults to \code{FALSE}.
#' @param matchvars Logical scalar specifying whether to compare the shared set of variable names in
#'   \code{comp} and \code{fossil} to the variable names in \code{struc} and pare them all down to the 
#'   set of shared variables. If \code{FALSE} and variable names differ then an error will be returned.
#'   Defaults to \code{FALSE}.
#' @param datastruc If multivariate data are used, this is a character string specifiying whether to 
#'   incorporate the missing data structure
#'   into dimorphism estimates (\code{"missing"}), whether to downsample to the missing data sample size but 
#'   keep all metric data for the comparative sample (\code{"complete"}), or to perform both types of
#'   resampling separately (\code{"both"}).  Ignored if only univariate data are provided or if all datasets
#'   are complete. Defaults to \code{NULL}, which reverts to \code{"missing"} if some datasets are incomplete.
#' @param templatevar A character object or integer value specifying the name or column number 
#'   of the variable in \code{fossil} and \code{comp} to be estimated using the template method.  Ignored 
#'   if template method is not used.  Defaults to \code{NULL}.
#' @return A list of class \code{SSDtest}.  Printing this object provides information about the datasets and 
#'   test performed, median values for resampled distributions for each sample for each method used, confidence 
#'   intervals for all estimators calculated, and calculated \emph{p}-values.  Plotting this object produces 
#'   histograms for resampled distributions for one or more estimation methods, or histograms for differences
#'   between samples (see \code{\link[dimorph]{plot.SSDtest}}).
#' @seealso \code{\link{pvals}}, \code{\link[dimorph]{plot.SSDtest}}
#' @references Gordon AD, Green DJ, Richmond BG. (2008) Strong postcranial size dimorphism in 
#'   \emph{Australopithecus afarensis}: Results from two new resampling methods for multivariate 
#'   data sets with missing data. \emph{American Journal of Physical Anthropology}. 135:311-328.
#'   (\href{https://doi.org/10.1002/ajpa.20745}{https://doi.org/10.1002/ajpa.20745}) 
#' @examples
#' # SSDtests using simulated fossils generated from real gorilla and human samples with some 
#' #   data removed
#' data(fauxil)
#' fauxil
#' 
#' ## Univariate examples
#' # First, some code that would generate errors
#' # SSDtest()
#' # SSDtest(fossil=list("Fauxil sp. 1"=fauxil[fauxil$Species=="Fauxil sp. 1", "FHSI"]))
#' # SSDtest(comp=list("G. gorilla"=apelimbart[apelimbart$Species=="Gorilla gorilla", "FHSI"]))
#' 
#' # Standard significance test with one fossil sample, sampling without replacement
#' test_faux_uni <- SSDtest(
#'   fossil=list("Fauxil sp. 1"=fauxil[fauxil$Species=="Fauxil sp. 1", "FHSI"]),
#'   comp=list("G. gorilla"=apelimbart[apelimbart$Species=="Gorilla gorilla", "FHSI"],
#'             "H. sapiens"=apelimbart[apelimbart$Species=="Homo sapiens", "FHSI"],
#'             "P. troglodytes"=apelimbart[apelimbart$Species=="Pan troglodytes", "FHSI"],
#'             "H. lar"=apelimbart[apelimbart$Species=="Hylobates lar", "FHSI"]),
#'   fossilsex=NULL,
#'   compsex=list("G. gorilla"=apelimbart[apelimbart$Species=="Gorilla gorilla", "Sex"],
#'                "H. sapiens"=apelimbart[apelimbart$Species=="Homo sapiens", "Sex"],
#'                "P. troglodytes"=apelimbart[apelimbart$Species=="Pan troglodytes", "Sex"],
#'                "H. lar"=apelimbart[apelimbart$Species=="Hylobates lar", "Sex"]),
#'   methsUni=c("SSD", "MMR", "BDI"),
#'   limit=1000,
#'   nResamp=100)
#' test_faux_uni
#' plot(test_faux_uni) # plots first method by default (SSD)
#' plot(test_faux_uni, est=2) # plots second method (MMR)
#' speciescolors <- c("Fauxil sp. 1"="#352A87", "Fauxil sp. 2"="#F9FB0E", "G. gorilla"="#EABA4B",
#'                    "H. sapiens"="#09A9C0", "P. troglodytes"="#78BE7C", "H. lar"="#0D77DA")
#' plot(test_faux_uni, est=2, groupcols=speciescolors) # change the colors
#' plot(test_faux_uni, type="diff", est=2) # plots differences between samples
#' plot(test_faux_uni, type="diff", est=2, diffs=c(1,2)) # plots diffs between first pair of samples
#' 
#' # Same as above, but variable name information is preserved
#' test_faux_uni2 <- SSDtest(
#'   fossil=list("Fauxil sp. 1"=fauxil[fauxil$Species=="Fauxil sp. 1", "FHSI", drop=FALSE]),
#'   comp=list("G. gorilla"=apelimbart[apelimbart$Species=="Gorilla gorilla", "FHSI", drop=FALSE],
#'             "H. sapiens"=apelimbart[apelimbart$Species=="Homo sapiens", "FHSI", drop=FALSE],
#'             "P. troglodytes"=apelimbart[apelimbart$Species=="Pan troglodytes","FHSI",drop=FALSE],
#'             "H. lar"=apelimbart[apelimbart$Species=="Hylobates lar", "FHSI", drop=FALSE]),
#'   fossilsex=NULL,
#'   compsex=list("G. gorilla"=apelimbart[apelimbart$Species=="Gorilla gorilla", "Sex"],
#'                "H. sapiens"=apelimbart[apelimbart$Species=="Homo sapiens", "Sex"],
#'                "P. troglodytes"=apelimbart[apelimbart$Species=="Pan troglodytes", "Sex"],
#'                "H. lar"=apelimbart[apelimbart$Species=="Hylobates lar", "Sex"]),
#'   methsUni=c("SSD", "MMR", "BDI"),
#'   limit=1000,
#'   nResamp=100)
#' test_faux_uni2
#' plot(test_faux_uni2, est=2, groupcols=speciescolors)
#' 
#' # Same as above except that null distributions are generated WITH replacement rather than WITHOUT
#' test_faux_uni3 <- SSDtest(
#'   fossil=list("Fauxil sp. 1"=fauxil[fauxil$Species=="Fauxil sp. 1", "FHSI", drop=FALSE]),
#'   comp=list("G. gorilla"=apelimbart[apelimbart$Species=="Gorilla gorilla", "FHSI", drop=FALSE],
#'             "H. sapiens"=apelimbart[apelimbart$Species=="Homo sapiens", "FHSI", drop=FALSE],
#'             "P. troglodytes"=apelimbart[apelimbart$Species=="Pan troglodytes","FHSI",drop=FALSE],
#'             "H. lar"=apelimbart[apelimbart$Species=="Hylobates lar", "FHSI", drop=FALSE]),
#'   fossilsex=NULL,
#'   compsex=list("G. gorilla"=apelimbart[apelimbart$Species=="Gorilla gorilla", "Sex"],
#'                "H. sapiens"=apelimbart[apelimbart$Species=="Homo sapiens", "Sex"],
#'                "P. troglodytes"=apelimbart[apelimbart$Species=="Pan troglodytes", "Sex"],
#'                "H. lar"=apelimbart[apelimbart$Species=="Hylobates lar", "Sex"]),
#'   methsUni=c("SSD", "MMR", "BDI"),
#'   limit=1000,
#'   nResamp=100,
#'   replace=TRUE)
#' test_faux_uni3
#' plot(test_faux_uni3, est=2, groupcols=speciescolors)
#' 
#' # Instead of standard test, each comparative and fossil sample is bootstrapped (sampled 
#' #   with replacement to that sample's full sample size) by setting 'fullsamplesboot' to
#' #   TRUE.  Produces a distribution of fossil estimates rather than a single point estimate, 
#' #   but note that comparative samples are resampled to their full sample size, not 
#' #   downsampled to the fossil sample size.
#' test_faux_uni4 <- SSDtest(
#'   fossil=list("Fauxil sp. 1"=fauxil[fauxil$Species=="Fauxil sp. 1", "FHSI", drop=FALSE]),
#'   comp=list("G. gorilla"=apelimbart[apelimbart$Species=="Gorilla gorilla", "FHSI", drop=FALSE],
#'             "H. sapiens"=apelimbart[apelimbart$Species=="Homo sapiens", "FHSI", drop=FALSE],
#'             "P. troglodytes"=apelimbart[apelimbart$Species=="Pan troglodytes","FHSI",drop=FALSE],
#'             "H. lar"=apelimbart[apelimbart$Species=="Hylobates lar", "FHSI", drop=FALSE]),
#'   fossilsex=NULL,
#'   compsex=list("G. gorilla"=apelimbart[apelimbart$Species=="Gorilla gorilla", "Sex"],
#'                "H. sapiens"=apelimbart[apelimbart$Species=="Homo sapiens", "Sex"],
#'                "P. troglodytes"=apelimbart[apelimbart$Species=="Pan troglodytes", "Sex"],
#'                "H. lar"=apelimbart[apelimbart$Species=="Hylobates lar", "Sex"]),
#'   methsUni=c("SSD", "MMR", "BDI"),
#'   limit=1000,
#'   nResamp=100,
#'   fullsamplesboot=TRUE)
#' test_faux_uni4
#' plot(test_faux_uni4, est=2, # plots the second estimation method (MMR)
#'      invert=1, # inverts the histogram for the first group in 'test_afar_uni4'
#'      groupcols=speciescolors)
#'
#' ## Multivariate examples
#' # GMM significance tests with a fossil sample for multiple estimators using both complete
#' #   and incomplete comparative datasets.  A single point estimate is generated for the fossil 
#' #   sample for each method.
#' SSDvars <- c("FHSI", "TPML", "TPMAP", "TPLAP", "HHMaj",
#'              "HHMin", "RHMaj", "RHMin", "RDAP", "RDML")
#' test_faux_multi1 <- SSDtest(
#'   fossil=list("Fauxil sp. 1"=fauxil[fauxil$Species=="Fauxil sp. 1", SSDvars]),
#'   comp=list("G. gorilla"=apelimbart[apelimbart$Species=="Gorilla gorilla", SSDvars],
#'             "H. sapiens"=apelimbart[apelimbart$Species=="Homo sapiens", SSDvars],
#'             "P. troglodytes"=apelimbart[apelimbart$Species=="Pan troglodytes", SSDvars],
#'             "H. lar"=apelimbart[apelimbart$Species=="Hylobates lar", SSDvars]),
#'   fossilsex=NULL,
#'   compsex=list("G. gorilla"=apelimbart[apelimbart$Species=="Gorilla gorilla", "Sex"],
#'                "H. sapiens"=apelimbart[apelimbart$Species=="Homo sapiens", "Sex"],
#'                "P. troglodytes"=apelimbart[apelimbart$Species=="Pan troglodytes", "Sex"],
#'                "H. lar"=apelimbart[apelimbart$Species=="Hylobates lar", "Sex"]),
#'   methsUni=c("MMR", "BDI"),
#'   methsMulti=c("GMM"),
#'   datastruc="both",
#'   nResamp=100)
#' test_faux_multi1
#' plot(test_faux_multi1, groupcols=speciescolors)
#' plot(test_faux_multi1, est=2, # plot the 2nd method combination: MMR, GMM, missing datastructure
#'      groupcols=speciescolors)
#' plot(test_faux_multi1, est=2,
#'      type="diff") # plot the diffs between samples for estimates using the 2nd method combination
#' 
#' # As above but with a different fossil sample.
#' test_faux_multi2 <- SSDtest(
#'   fossil=list("Fauxil sp. 2"=fauxil[fauxil$Species=="Fauxil sp. 2", SSDvars]),
#'   comp=list("G. gorilla"=apelimbart[apelimbart$Species=="Gorilla gorilla", SSDvars],
#'             "H. sapiens"=apelimbart[apelimbart$Species=="Homo sapiens", SSDvars],
#'             "P. troglodytes"=apelimbart[apelimbart$Species=="Pan troglodytes", SSDvars],
#'             "H. lar"=apelimbart[apelimbart$Species=="Hylobates lar", SSDvars]),
#'   fossilsex=NULL,
#'   compsex=list("G. gorilla"=apelimbart[apelimbart$Species=="Gorilla gorilla", "Sex"],
#'                "H. sapiens"=apelimbart[apelimbart$Species=="Homo sapiens", "Sex"],
#'                "P. troglodytes"=apelimbart[apelimbart$Species=="Pan troglodytes", "Sex"],
#'                "H. lar"=apelimbart[apelimbart$Species=="Hylobates lar", "Sex"]),
#'   methsUni=c("MMR", "BDI"),
#'   methsMulti=c("GMM"),
#'   datastruc="both",
#'   nResamp=100)
#' test_faux_multi2
#' plot(test_faux_multi2, groupcols=speciescolors)
#' plot(test_faux_multi2, est=2, groupcols=speciescolors)
#' plot(test_faux_multi2, est=2, type="diff")
#' 
#' # Now the same data using a much more conservative approach by setting 'rebootstrap' to TRUE.
#' test_faux_multi3 <- SSDtest(
#'   fossil=list("Fauxil sp. 2"=fauxil[fauxil$Species=="Fauxil sp. 2", SSDvars]),
#'   comp=list("G. gorilla"=apelimbart[apelimbart$Species=="Gorilla gorilla", SSDvars],
#'             "H. sapiens"=apelimbart[apelimbart$Species=="Homo sapiens", SSDvars],
#'             "P. troglodytes"=apelimbart[apelimbart$Species=="Pan troglodytes", SSDvars],
#'             "H. lar"=apelimbart[apelimbart$Species=="Hylobates lar", SSDvars]),
#'   fossilsex=NULL,
#'   compsex=list("G. gorilla"=apelimbart[apelimbart$Species=="Gorilla gorilla", "Sex"],
#'                "H. sapiens"=apelimbart[apelimbart$Species=="Homo sapiens", "Sex"],
#'                "P. troglodytes"=apelimbart[apelimbart$Species=="Pan troglodytes", "Sex"],
#'                "H. lar"=apelimbart[apelimbart$Species=="Hylobates lar", "Sex"]),
#'   methsUni=c("MMR", "BDI"),
#'   methsMulti=c("GMM"),
#'   datastruc="both",
#'   nResamp=100,
#'   rebootstrap=TRUE)
#' test_faux_multi3
#' plot(test_faux_multi3, est=2, invert=1, groupcols=speciescolors)
#' @export
SSDtest <- function(fossil=NULL,
                    comp=NULL,
					fossilsex=NULL,
					compsex=NULL,
					#H0=NULL,
					methsUni=c("MMR", "BDI"),
					methsMulti="GMM",
					replace=F,
					rebootstrap=F,
					fullsamplesboot=F, # if all complete and this is true, bootstraps each sample
					nResamp=1000,
					exactcomp=T,
					exactfossil=T,
					limit=50000,
					center="geomean",
					sex.female=1,
					na.rm=T,
					ncorrection=F,
					matchvars=F,
					datastruc=NULL,
					templatevar=NULL) {
  # fossil or comp can be NULL, but there must be at least two samples to run the
  #   test (but it could be a comparison just between fossils or just between
  #   comparative data).
  #
  # I think I may want to remove the H0 option.  It really only makes sense to use an H0
  #   in two scenrarios: the first is to test against monomorphism for a sample with no
  #   missing data, in which case 'bootdimorph' can be used to generate a confidence interval
  #   and the limits can be checked against zero.  The second is to check against the level of
  #   dimorphism observed in a fossil sample, in which case you need to incorporate the missing
  #   data structure of that fossil sample to make sure that dimorphism is estimated correctly.
  #   So including an H0 option is not useful here, and it could lead to incorrectly testing
  #   the level of dimorphism in a complete sample against a fossil point estimate derived from
  #   a missing data structure.
  #
  # If list, names of list elements can be included to name samples.
  # If list, order samples are provided will be order one-sided tests will be generated (most dimorphic first).
  # Maybe also code so that x can be an object previously produced by SSDtest,
  #   and is just being rerun to calculate p values for different SSD estimate methods.
  #   In that case, everything except methsUni, methsMulti, and templatevar would
  #   be ignored.
  # compsex can either be a list of same length as x or a vector of length equal to the number of 
  #   rows in x - contains sex information for each sample in x.
  # In results, order samples by median SSD estimate (first specified method in methsUni).
  # Note: can't use template method with vector structure
  methsUni <- match.arg(methsUni, choices=c("SSD", "MMR", "BDI", "MoM", "FMA", "BFM",
                                            "ERM", "CV", "CVsex", "sdlog", "sdlogsex"),
						several.ok=T)
  methsMulti <- match.arg(methsMulti, choices=c("GMM", "GMsize", "TM"), several.ok=T)
  center <- match.arg(center, choices=c("geomean", "mean"), several.ok=T)
  datastruc <- match.arg(datastruc, choices=c("missing", "complete", "both"))
  fossilNull <- TRUE # Flag for running standard randomization test
  if (!is.null(fossil)) {
    fossilNull <- FALSE
    if (!inherits(fossil, c("list", "data.frame", "matrix"))) stop("'fossil' must be a data frame or matrix holding metric data for a sample, or a list of such objects.")
    #if (!class(fossil) %in% c("list", "data.frame", "matrix")) stop("'fossil' must be a data frame or matrix holding metric data for a sample, or a list of such objects.")
    if (!inherits(fossil, "list")) fossil <- list(fossil)
    #if (class(fossil)!="list") fossil <- list(fossil)
  }
  if (!is.null(comp)) {
    if (!inherits(comp, c("list", "data.frame", "matrix"))) stop("'comp' must be a data frame or matrix holding metric data for a sample, or a list of such objects.")
    #if (!class(comp) %in% c("list", "data.frame", "matrix")) stop("'comp' must be a data frame or matrix holding metric data for a sample, or a list of such objects.")
    if (!inherits(comp, "list")) comp <- list(comp)
    #if (class(comp)!="list") comp <- list(comp)
  }
  #if (!class(x) %in% c("list", "data.frame", "matrix")) stop("'x' must be a data frame or matrix holding metric data for a sample, or a list of such objects.")
  #if (class(x)!="list") x <- list(x)
  nfoss <- 0
  ncomp <- 0
  if (!is.null(fossil)) nfoss <- length(fossil)
  if (!is.null(comp)) ncomp <- length(comp)
  ndat <- nfoss + ncomp
  if (ndat < 2) stop("'SSDtest' requires at least two samples.")
  if (is.null(fossil)) warning("No fossil data have been designated - proceeding assuming all data are comparative.")
  if (is.null(comp)) warning("No comparative data have been designated - proceeding assuming all data are fossil.")
  if (!is.null(fossil)) if (is.null(names(fossil))) names(fossil) <- paste0("fossil_sample_", 1:nfoss)
  if (!is.null(comp)) if (is.null(names(comp))) names(comp) <- paste0("comparative_sample_", 1:ncomp)
  #if (is.null(names(x))) names(x) <- paste0("sample_", 1:ndat)
  # check each element of fossil and convert to dataframe if vector or matrix
  if (!is.null(fossil)) for (i in 1:nfoss) {
    uni <- NA
    {if ("data.frame" %in% class(fossil[[i]])| "matrix" %in% class(fossil[[i]])) {
	  if ("matrix" %in% class(fossil[[i]])) fossil[[i]] <- data.frame(fossil[[i]])
      {if (ncol(fossil[[i]])==1) uni <- T
      else if (ncol(fossil[[i]])>1) uni <- F}  
    }
    else if ("numeric" %in% class(fossil[[i]])| "integer" %in% class(fossil[[i]])) {
	  uni <- T
	  fossil[[i]] <- data.frame(fossil[[i]])
	  vars <- "VAR"
	  colnames(fossil[[i]]) <- vars
	  methsMulti <- NULL
	  templatevar <- NULL
    }
	}
    if (is.na(uni)) stop("paste0('", names(fossil)[[i]], "' must be a numeric vector, matrix, or dataframe.")
	# remove empty cases
	fossilkeep <- apply(fossil[[i]], 1, function (y) return(sum(is.na(y)))) < ncol(fossil[[i]])
    if (prod(fossilkeep)==0) {
      fossil[[i]] <- fossil[[i]][fossilkeep,,drop=F]
	  if(!is.null(fossilsex[[i]])) fossilsex[[i]] <- fossilsex[[i]][fossilkeep]
    }
  }
  # check each element of comp and convert to dataframe if vector or matrix
  # also, check that each comp sample has no missing data
  if (!is.null(comp)) for (i in 1:ncomp) {
    uni <- NA
    {if ("data.frame" %in% class(comp[[i]])| "matrix" %in% class(comp[[i]])) {
	  if ("matrix" %in% class(comp[[i]])) comp[[i]] <- data.frame(comp[[i]])
      {if (ncol(comp[[i]])==1) uni <- T
      else if (ncol(comp[[i]])>1) uni <- F}  
    }
    else if ("numeric" %in% class(comp[[i]])| "integer" %in% class(comp[[i]])) {
	  uni <- T
	  comp[[i]] <- data.frame(comp[[i]])
	  vars <- "VAR"
	  colnames(comp[[i]]) <- vars
	  methsMulti <- NULL
	  templatevar <- NULL
    }
	}
    if (is.na(uni)) stop("paste0('", names(comp)[[i]], "' must be a numeric vector, matrix, or dataframe.")
	# check for incomplete cases
	compcomplete <- stats::complete.cases(comp[[i]])
    if (prod(compcomplete)==0) {
      {if (na.rm) {
	    comp[[i]] <- comp[[i]][compcomplete,,drop=F]
	    if(!is.null(compsex[[i]])) compsex[[i]] <- compsex[[i]][compcomplete]
		warning(paste0(names(comp)[i], " in 'comp' has missing data: incomplete specimens removed from analysis."))
	  }
	  else stop("At least one sample in 'comp' has missing data and 'na.rm'=FALSE.")}
	}
  }
  # check fossilsex information
  if (!is.null(fossilsex)) {
    if (length(fossilsex)!=nfoss) stop("'fossilsex' must either be NULL if no sex information is present for fossil samples, or a list\ncontaining either NULL or sex information for each sample (see 'dimorph'\nfor more information on formatting sex information).")
	if (is.null(names(fossilsex))) names(fossilsex) <- names(fossil)
	if (!identical(names(fossil), names(fossilsex))) stop("The names of list elements for 'fossil' and 'fossilsex' must match. Alternatively,\nif the names of 'fossilsex' elements are unspecified then they will assume\nto be provided in the same order as the elements of 'fossil'.")
    for (i in 1:nfoss) {
	  if (!is.null(fossilsex[[i]])) if(nrow(fossil[[i]])!=length(fossilsex[[i]])) stop("When sex information is provided, the number of individuals in the sex\nvector in 'fossilsex' must match the number of individuals in the corresponding\ndata set in 'fossil'.")
	}
  }
  # check compsex information
  if (!is.null(compsex)) {
    if (length(compsex)!=ncomp) stop("'compsex' must either be NULL if no sex information is present, or a list\ncontaining either NULL or sex information for each sample (see 'dimorph'\nfor more information on formatting sex information).")
	if (is.null(names(compsex))) names(compsex) <- names(comp)
	if (!identical(names(comp), names(compsex))) stop("The names of list elements for 'comp' and 'compsex' must match. Alternatively,\nif the names of 'compsex' elements are unspecified then they will assume\nto be provided in the same 	order as the elements of 'comp'.")
    for (i in 1:ncomp) {
	  if (!is.null(compsex[[i]])) if(nrow(comp[[i]])!=length(compsex[[i]])) stop("When sex information is provided, the number of individuals in the sex\nvector in 'compsex' must match the number of individuals in the corresponding\ndata set in 'comp'.")
	}
  }
  #################################
  # STANDARD RANDOMIZATION IF NO FOSSIL SAMPLE
  #################################
  if (fossilNull) { # run standard randomization and return object
    out <- NULL
    # add warning that a standard randomization will be performed and some arguments may be ignored
    # add option to constrain sampling to same sex? i.e., randomize within males and randomize within females.
    # 'estimates' will be a list where each element is not a sample, but rather a pairwise comparison 
    #    between samples
    # for each pair, calculate number of exact combinations if exact is selected and compare to 'limit',
    #    running Monte Carlo if necessary
    # for each pair, record the smaller sample name (and choose the lower factor number if they are tied
    #    for sample size) and then record the adresses chosen for the smaller sample from a combined vector
    #    of both samples for each iteration, where the first iteration is the actual sample.  This will 
    #    necessitate recording different sampling addresses for every pair.
    # 'estimate' data frame for each pair will contain an estimate of the difference for that pair rather than 
    #    an estimate of the value for a sample.  A column should be added to denot which set of addresses was used
    #    for a particular difference, with the first set of addresses always being the actual difference.
  
  
    return(out)
  }
  #################################
  # END OF STANDARD RANDOMIZATION IF NO FOSSIL SAMPLE
  #################################
  # get structure, but do differently if all data are complete and
  #  it's just a test using full sample sizes for all samples,
  #  in which case I should be using bootdimorph to get the
  #  resampled values.
  {if (!is.null(fossil)) {
    struct <- dimorph::getstructure(fossil, forcematrix=T)
    {if (inherits(struct, "data.frame")) {
    #{if (class(struct)=="data.frame") {
	  vars <- colnames(struct)
	  nfossspec <- nrow(struct)
    }
	else {
      vars <- names(struct)
      nfossspec <- max(struct)
	  if (!is.null(methsMulti)) {
	    if ("TM" %in% methsMulti) {
	      warning("The template method cannot be used when two or more fossil samples have missing data.")
		  methsMulti <- methsMulti[!methsMulti %in% "TM"]
		  if (length(methsMulti)==0) stop("No user-specified multivariate methods remain after removal of template method.")
	    }
	  }
    }}
	for (i in 1:nfoss) {
	  fossil[[i]] <- fossil[[i]][,vars,drop=F]
	  # remove rows with no data
	  fkeep <- apply(fossil[[i]], 1, function (x) return(sum(is.na(x)))) < length(vars)
	  fossil[[i]] <- fossil[[i]][fkeep,,drop=F]
	  if (!is.null(fossilsex[[i]])) fossilsex[[i]] <- fossilsex[[i]][fkeep]
	}
	if (!is.null(comp)) {
	  for (i in 1:ncomp) {
	    if (prod(vars %in% colnames(comp[[i]]))==0) stop("Not all variables shared by fossil sample(s) present in all comparative samples.")
	    comp[[i]] <- comp[[i]][,vars,drop=F]
	  }
      if (!prod(unlist(lapply(comp, nrow)) > nfossspec)) {
	    if (replace) warning("At least one comparative sample has fewer specimens than at least one fossil sample.")
		else stop("At least one comparative sample has fewer specimens than at least one fossil sample.")
	  }  
	}	
  }
  else {
    struct <- dimorph::getstructure(comp)
	vars <- colnames(struct)
	for (i in 1:ncomp) comp[[i]] <- comp[[i]][,vars,drop=F]
  }}
  complete <- unlist(lapply(c(fossil,comp), function(x) as.logical(prod(stats::complete.cases(x)))))
  {if (prod(complete)==1 & fullsamplesboot) { # if all data sets are complete & fullsamplesboot arg is TRUE
    # use bootdimorph here unless otherwise specified
	xSSDres <- list(NULL)
	i <- 0
	if (!is.null(fossil)) for (i in 1:length(fossil)) {
	  strc <- NULL
	  datast <- "complete"
	  # if (!complete[i]) datast <- "missing" # shouldn't ever get triggered
	  xSSDres[[i]] <- dimorph::bootdimorph(fossil[[i]],
	                              sex=fossilsex[[i]],
								  methsUni=methsUni,
								  methsMulti=methsMulti,
								  nResamp=nResamp,
								  exact=exactcomp,
								  limit=limit,
								  center=center,
								  sex.female=sex.female,
								  na.rm=na.rm,
								  ncorrection=ncorrection,
								  struc=strc,
								  #datastruc=NULL,
								  datastruc=datast,
								  templatevar=templatevar,
								  alternative="two.sided",
								  conf.level=0.95)
	  attr(xSSDres[[i]], "fullsamplesboot") <- T 
	  names(xSSDres)[i] <- names(fossil)[i]
	}
	j <- i
	if (!is.null(comp)) for (i in 1:length(comp)) {
	  strc <- NULL
	  datast <- "complete"
	  xSSDres[[i+j]] <- dimorph::bootdimorph(comp[[i]],
	                              sex=compsex[[i]],
								  methsUni=methsUni,
								  methsMulti=methsMulti,
								  nResamp=nResamp,
								  exact=exactcomp,
								  limit=limit,
								  center=center,
								  sex.female=sex.female,
								  na.rm=na.rm,
								  ncorrection=ncorrection,
								  struc=strc,
								  #datastruc=NULL,
								  datastruc=datast,
								  templatevar=templatevar,
								  alternative="two.sided",
								  conf.level=0.95)
	  attr(xSSDres[[i+j]], "fullsamplesboot") <- T 
	  names(xSSDres)[i+j] <- names(comp)[i]
	}
  } # end if all data sets are complete & fullsamplesboot arg is TRUE
  else { # if at least one data set has missing data or if fullsamplesboot is FALSE
    # if only one fossil sample then we need to estimate dimorphism as a point estimate for that sample.
    xAds <- list(NULL)
	{if (nfoss==1) { # only one or zero missing data set and if one, it's equal to structure
	  # For single fossil sample, builds a dimorphAds object assuming the dataset is complete,
	  #   then replaces the relevant parts of the object to reflect that it's actually incomplete
	  #   if it is incomplete.  Generates a single set of addresses using all data for a point estimate.
	  tmpdat <- fossil[[1]]
	  tmpdat[is.na(tmpdat)] <- 1
	  xAds[[1]] <- dimorph::getsampleaddresses(comparative=tmpdat,
	                                struc=tmpdat,
									compsex=fossilsex[[1]],
									nResamp=1,
									exact=F,
									matchvars=matchvars,
									replace=F)
									#replace=replace)
	  xAds[[1]]$comparative <- fossil[[1]]
	  xAds[[1]]$struc <- struct
	  xAds[[1]]$addresses[,1] <- 1:nrow(xAds[[1]]$addresses)
	  xAds[[1]]$exact <- T
	  # Now it gets the resampled addresses for the comparative data
	  for (i in 1:length(comp)) {
	    xAds[[i+1]] <- dimorph::getsampleaddresses(comparative=comp[[i]],
	                                struc=struct,
									compsex=compsex[[i]],
									nResamp=nResamp,
									exact=exactcomp,
									limit=limit,
									matchvars=matchvars,
									#replace=F)
									replace=replace)
	  }
	} # end only one fossil
	else { # other than one fossil
	  if (!is.null(fossil)) for (i in 1:length(fossil)) {
	    xAds[[i]] <- dimorph::getsampleaddresses(comparative=fossil[[i]],
	                                struc=struct,
									compsex=fossilsex[[i]],
									nResamp=nResamp,
									exact=exactfossil,
									limit=limit,
									matchvars=matchvars,
									replace=F)
									#replace=replace)
      }
	  if (!is.null(comp)) for (i in 1:length(comp)) {
	    xAds[[i+nfoss]] <- dimorph::getsampleaddresses(comparative=comp[[i]],
	                                struc=struct,
									compsex=compsex[[i]],
									nResamp=nResamp,
									exact=exactcomp,
									limit=limit,
									matchvars=matchvars,
									#replace=F)
									replace=replace)
      }
	} # end other than one fossil
	} # end if else
	names(xAds) <- c(names(fossil), names(comp))
	if (rebootstrap) {
    # add procedure to run rebootstrap for every sample; how do 
	#   I want to caclulate or specify the nResamp for this function,
	#   given that some samples may have used exact resampling and some might not?
	# In any event, this will return an updated version of xAds.
	  for (i in 1:ndat) {
	    nR <- NA
		if (xAds[[i]]$nResamp==1) nR <- nResamp
	    xAds[[i]] <- dimorph:::bootsampleaddresses(xAds[[i]], nResamp=nR)
      }
	}
	# Now run resampleSSD
	xSSDres <- list(NULL)
	if (!is.null(fossil)) for (i in 1:nfoss) {
	  # figure out when to include sex-dependent methods and when not to.
	  methsUnitmp <- methsUni
	  if (is.null(fossilsex[[i]])) methsUnitmp <- methsUnitmp[!methsUnitmp %in% c("SSD", "CVsex", "sdlogsex")]
	  xSSDres[[i]] <- dimorph::resampleSSD(xAds[[i]],
	                            #datastruc="missing",
	                            datastruc=datastruc,
								methsMulti=methsMulti,
								methsUni=methsUnitmp,
								sex.female=sex.female,
								center=center,
								templatevar=templatevar,
								na.rm=T,
								ncorrection=ncorrection,
								details=F)
	  names(xSSDres)[i] <- names(fossil)[i]
    }
	if (!is.null(comp)) for (i in 1:ncomp) {
	  # figure out when to include sex-dependent methods and when not to.
	  methsUnitmp <- methsUni
	  if (is.null(compsex[[i]])) methsUnitmp <- methsUnitmp[!methsUnitmp %in% c("SSD", "CVsex", "sdlogsex")]
	  xSSDres[[i+nfoss]] <- dimorph::resampleSSD(xAds[[i+nfoss]],
	                            #datastruc="missing",
	                            datastruc=datastruc,
								methsMulti=methsMulti,
								methsUni=methsUnitmp,
								sex.female=sex.female,
								center=center,
								templatevar=templatevar,
								na.rm=T,
								ncorrection=ncorrection,
								details=F)
	  names(xSSDres)[i+nfoss] <- names(comp)[i]
    }
  } # end if at least one data set has missing data or if fullsamplesboot is FALSE
  } # end if else
  # - get full set of unique analysis types (combination of methUni, methMulti, datastructure, and center)
  #   used across all samples.
  datnames <- c(names(fossil), names(comp))
  multi <- F
  if (ncol(xSSDres[[1]]$sampleADS$comparative) > 1) multi <- T
  {if (multi) {
    methcombos <- levels(droplevels(xSSDres[[1]]$estimates$methodUni:xSSDres[[1]]$estimates$methodMulti:xSSDres[[1]]$estimates$center:xSSDres[[1]]$estimates$datastructure))
	if (ndat > 1) for (i in 2:ndat) {
	  methcombostmp <- levels(droplevels(xSSDres[[i]]$estimates$methodUni:xSSDres[[i]]$estimates$methodMulti:xSSDres[[i]]$estimates$center:xSSDres[[i]]$estimates$datastructure))
	  methcombos <- unique(c(methcombostmp, methcombos))
	}
	methcombos <- data.frame(do.call(rbind, strsplit(methcombos, ":")))
	colnames(methcombos) <- c("methodUni", "methodMulti", "center", "datastructure") 
    # - figure out which samples have estimates for each analysis type (e.g., no SSD for samples without sex info).
    # - for each analysis type, calculate median and mean estimate for each sample.
    for (i in 1:ndat) {
      methcombos[,paste0("median.", datnames[i])] <- NA
      methcombos[,paste0("mean.", datnames[i])] <- NA
      for (j in 1:nrow(methcombos)) {
	    tmp <- xSSDres[[i]]$estimates$estimate[xSSDres[[i]]$estimates$methodUni==methcombos$methodUni[j] &
	                                         xSSDres[[i]]$estimates$methodMulti==methcombos$methodMulti[j] &
											 xSSDres[[i]]$estimates$center==methcombos$center[j] &
											 xSSDres[[i]]$estimates$datastructure==methcombos$datastructure[j]]
	    if (length(tmp) > 0) {
	      methcombos[j,paste0("median.", datnames[i])] <- median(tmp, na.rm=T)
	      methcombos[j,paste0("mean.", datnames[i])] <- mean(tmp, na.rm=T)
	    }
	  }
    }
    # - for each analysis type, within each possible pair of samples, calculate all possible differences
    #   in the estimate between one value for one sample and one value for the other sample, using the order in
    #   which samples were provided in the list 'x'.  Store these in a list of lists: first level of list is the unique 
    #   analysis type; second level of list is vector of differences for a particular pair of samples. If H0
    #   is provided, then rather than calculating p for pairwise differences, calculate p for difference of 
    #   each sample from the H0 value for that univariate method.
    {if (ndat > 1) {
      pair <- data.frame(sample.greater=rep(NA, choose(ndat,2)), sample.lesser=rep(NA, choose(ndat,2)))
	  counter <- 0
      for (i in 1:(ndat-1)) {
	    for (j in (i+1):ndat) {
	      counter <- counter+1
		  pair$sample.greater[counter] <- datnames[i]
		  pair$sample.lesser[counter] <- datnames[j]
	    }
	  }
	  xDiffs <- list(NULL)
	  for (i in 1:nrow(methcombos)) {
	    xDiffs[[i]] <- list(NULL)
	    names(xDiffs)[i] <- paste0(methcombos$methodUni[i], ":",
	                           methcombos$methodMulti[i], ":",
	                           methcombos$center[i], ":",
	                           methcombos$datastructure[i])
	    #{if (is.null(H0)) {
	      for (j in 1:nrow(pair)) {
		    grt <- xSSDres[[pair$sample.greater[j]]]$estimates$estimate[xSSDres[[pair$sample.greater[j]]]$estimates$methodUni==methcombos$methodUni[i] &
	                        xSSDres[[pair$sample.greater[j]]]$estimates$methodMulti==methcombos$methodMulti[i] &
							xSSDres[[pair$sample.greater[j]]]$estimates$center==methcombos$center[i] &
							xSSDres[[pair$sample.greater[j]]]$estimates$datastructure==methcombos$datastructure[i]]
		    les <- xSSDres[[pair$sample.lesser[j]]]$estimates$estimate[xSSDres[[pair$sample.lesser[j]]]$estimates$methodUni==methcombos$methodUni[i] &
	                        xSSDres[[pair$sample.lesser[j]]]$estimates$methodMulti==methcombos$methodMulti[i] &
							xSSDres[[pair$sample.lesser[j]]]$estimates$center==methcombos$center[i] &
							xSSDres[[pair$sample.lesser[j]]]$estimates$datastructure==methcombos$datastructure[i]]
		    xDiffs[[i]][[j]] <- rep(grt, each=length(les)) - rep(les, length(grt))
		    names(xDiffs[[i]])[j] <- paste0(pair$sample.greater[j], " - ", pair$sample.lesser[j])
	      }
		  for (j in 1:nrow(pair)) {
		    xDiffs[[i]][[j+nrow(pair)]] <- -1 * xDiffs[[i]][[j]]
			names(xDiffs[[i]])[j+nrow(pair)] <- paste0(pair$sample.lesser[j], " - ", pair$sample.greater[j])
		  }
	    #} # end is.null(H0)
	    #else { # H0 is not null
	    #  for (j in 1:ndat) {
		#    xDiffs[[i]][[j]] <- xSSDres[[j]]$estimates$estimate[xSSDres[[j]]$estimates$methodUni==methcombos$methodUni[i] &
	    #                                     xSSDres[[j]]$estimates$methodMulti==methcombos$methodMulti[i] &
		#									 xSSDres[[j]]$estimates$center==methcombos$center[i] &
		#									 xSSDres[[j]]$estimates$datastructure==methcombos$datastructure[i]] - H0[methcombos$methodUni[i]] # H0 with name of methodUni
		#    names(xDiffs[[i]])[j] <- names(x)[j]
		#  }
	    #} # end H0 is not null
	    #} # end if else
      }
    } # end ndat > 1
    else { # ndat==1: For one sample test, instead use H0 value to calculate difference from.
	  # This should never trigger now, since ndat must be 2 or higher
	  #xDiffs <- list(NULL)
	  #for (i in 1:nrow(methcombos)) {
	  #  xDiffs[[i]] <- list(NULL)
	  #  names(xDiffs)[i] <- paste0(methcombos$methodUni[i], ":",
	  #                         methcombos$methodMulti[i], ":",
	  #                         methcombos$center[i], ":",
	  #                         methcombos$datastructure[i])
	  #  xDiffs[[i]][[1]] <- xSSDres[[1]]$estimates$estimate[xSSDres[[1]]$estimates$methodUni==methcombos$methodUni[i] &
	  #                          xSSDres[[1]]$estimates$methodMulti==methcombos$methodMulti[i] &
	  #							xSSDres[[1]]$estimates$center==methcombos$center[i] &
	  #							xSSDres[[1]]$estimates$datastructure==methcombos$datastructure[i]] - H0[methcombos$methodUni[i]] # H0 with name of methodUni
	  #	names(xDiffs[[i]])[1] <- datnames[1]
      #} 
    } # end ndat==1
    } # end if else
  } # end if multi
  else { # univariate analysis
    methcombos <- levels(droplevels(xSSDres[[1]]$estimates$methodUni:xSSDres[[1]]$estimates$center))
	if (ndat > 1) for (i in 2:ndat) {
	  methcombostmp <- levels(droplevels(xSSDres[[i]]$estimates$methodUni:xSSDres[[i]]$estimates$center))
	  methcombos <- unique(c(methcombostmp, methcombos))
	}
	methcombos <- data.frame(do.call(rbind, strsplit(methcombos, ":")))
	colnames(methcombos) <- c("methodUni", "center") 
    # - figure out which samples have estimates for each analysis type (e.g., no SSD for samples without sex info).
    # - for each analysis type, calculate median and mean estimate for each sample.
    for (i in 1:ndat) {
      methcombos[,paste0("median.", datnames[i])] <- NA
      methcombos[,paste0("mean.", datnames[i])] <- NA
      for (j in 1:nrow(methcombos)) {
	    tmp <- xSSDres[[i]]$estimates$estimate[xSSDres[[i]]$estimates$methodUni==methcombos$methodUni[j] &
											 xSSDres[[i]]$estimates$center==methcombos$center[j]]
	    if (length(tmp) > 0) {
	      methcombos[j,paste0("median.", datnames[i])] <- median(tmp, na.rm=T)
	      methcombos[j,paste0("mean.", datnames[i])] <- mean(tmp, na.rm=T)
	    }
	  }
    }
    # - for each analysis type, within each possible pair of samples, calculate all possible differences
    #   in the estimate between one value for one sample and one value for the other sample, using the order in
    #   which samples were provided in the list 'x'.  Store these in a list of lists: first level of list is the unique 
    #   analysis type; second level of list is vector of differences for a particular pair of samples. If H0
    #   is provided, then rather than calculating p for pairwise differences, calculate p for difference of 
    #   each sample from the H0 value for that univariate method.
    {if (ndat > 1) {
      pair <- data.frame(sample.greater=rep(NA, choose(ndat,2)), sample.lesser=rep(NA, choose(ndat,2)))
      #pair <- data.frame(sample.greater=rep(NA, choose(ndat,2)*2), sample.lesser=rep(NA, choose(ndat,2)*2))
	  counter <- 0
      for (i in 1:(ndat-1)) {
	    for (j in (i+1):ndat) {
	      counter <- counter+1
		  pair$sample.greater[counter] <- datnames[i]
		  pair$sample.lesser[counter] <- datnames[j]
	    }
	  }
	  xDiffs <- list(NULL)
	  for (i in 1:nrow(methcombos)) {
	    xDiffs[[i]] <- list(NULL)
	    names(xDiffs)[i] <- paste0(methcombos$methodUni[i], ":",
	                           methcombos$center[i])
	    #{if (is.null(H0)) {
	      for (j in 1:nrow(pair)) {
		    grt <- xSSDres[[pair$sample.greater[j]]]$estimates$estimate[xSSDres[[pair$sample.greater[j]]]$estimates$methodUni==methcombos$methodUni[i] &
						xSSDres[[pair$sample.greater[j]]]$estimates$center==methcombos$center[i]]
		    les <- xSSDres[[pair$sample.lesser[j]]]$estimates$estimate[xSSDres[[pair$sample.lesser[j]]]$estimates$methodUni==methcombos$methodUni[i] &
						xSSDres[[pair$sample.lesser[j]]]$estimates$center==methcombos$center[i]]
		    xDiffs[[i]][[j]] <- rep(grt, each=length(les)) - rep(les, length(grt))
		    names(xDiffs[[i]])[j] <- paste0(pair$sample.greater[j], " - ", pair$sample.lesser[j])
	      }
		  for (j in 1:nrow(pair)) {
		    xDiffs[[i]][[j+nrow(pair)]] <- -1 * xDiffs[[i]][[j]]
			names(xDiffs[[i]])[j+nrow(pair)] <- paste0(pair$sample.lesser[j], " - ", pair$sample.greater[j])
		  }
	    #} # end is.null(H0)
	    #else { # H0 is not null
	    #  for (j in 1:ndat) {
		#    xDiffs[[i]][[j]] <- xSSDres[[j]]$estimates$estimate[xSSDres[[j]]$estimates$methodUni==methcombos$methodUni[i] &
		#									 xSSDres[[j]]$estimates$center==methcombos$center[i]] - H0[methcombos$methodUni[i]] # H0 with name of methodUni
		#    names(xDiffs[[i]])[j] <- names(x)[j]
		#  }
	    #} # end H0 is not null
	    #} # end if else
      }
    } # end ndat > 1
    else { # ndat==1: For one sample test, instead use H0 value to calculate difference from.
	  # This should never trigger now, since ndat must be 2 or higher
	  #xDiffs <- list(NULL)
	  #for (i in 1:nrow(methcombos)) {
	  #  xDiffs[[i]] <- list(NULL)
	  #  names(xDiffs)[i] <- paste0(methcombos$methodUni[i], ":", methcombos$center[i])
	  #  xDiffs[[i]][[1]] <- xSSDres[[1]]$estimates$estimate[xSSDres[[1]]$estimates$methodUni==methcombos$methodUni[i] &
	  #							xSSDres[[1]]$estimates$center==methcombos$center[i]] - H0[methcombos$methodUni[i]] # H0 with name of methodUni
	  #	names(xDiffs[[i]])[1] <- datnames[1]
      #} 
    } # end ndat==1
    } # end if else
  } # end univariate
  } # end if else
  # - For one-sided test, calculate proportion of differences less than or equal to zero. For two-sided
  #   test, also calculate proportion of differences that are that far from the mean or farther in the other
  #   direction (i.e., proportion of differences that are two times the mean or greater).
  pone <- function(diffs) {
    diffs <- diffs[!is.na(diffs)]
	if (length(diffs)==0) return(NA)
	return(sum(diffs <= 0)/length(diffs))
  }
  ptwo <- function(diffs) {
    diffs <- diffs[!is.na(diffs)]
	if (length(diffs)==0) return(NA)
	p <- (sum(diffs <= 0) + sum(diffs >= 2*mean(diffs)))/length(diffs)
	if (p > 1) {
	  p <- (sum(-diffs <= 0) + sum(-diffs >= 2*mean(-diffs)))/length(diffs)
	}
	return(p)
  }
  p.onesided <- do.call(cbind, lapply(xDiffs, function(y) do.call(rbind, lapply(y, pone))))
  colnames(p.onesided) <- names(xDiffs)
  p.twosided <- do.call(cbind, lapply(xDiffs, function(y) do.call(rbind, lapply(y, ptwo))))
  colnames(p.twosided) <- names(xDiffs)
  p.twosided <- p.twosided[1:(nrow(p.twosided)/2),] # same values for A-B and B-A
  p <- list(p.onesided=p.onesided, p.twosided=p.twosided)
  wrntxt <- NULL
  for (i in 1:length(xDiffs)) {
    for (j in 1:length(xDiffs[[i]])) {
      if (sum(is.na(xDiffs[[i]][[j]])) > 0) wrntxt <- paste0(wrntxt, "  ", names(xDiffs)[i], ": ",  names(xDiffs[[i]])[j], "\n")
	}
  }
  if (!is.null(wrntxt)) {
    wrntxt <- paste0("The following comparisons contain NAs that were\ndropped in the calculation of p-values:\n", wrntxt)
	warning(wrntxt)
  }
  out <- list(estimates=xSSDres,
			  #differences=xDiffs, # do I really need to save these? These can be recaluclated when needed.
			  methcombos=methcombos,
			  #H0=H0,
			  pvalues=p)
  class(out) <- "SSDtest"
  attr(out, "warntext") <- wrntxt
  return(out)
}


