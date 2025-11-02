#' @noRd
get.objects <- function() {
  # returns table of objects in memory and their sizes in bytes
  object.sizes <- NULL
  object.sizes <- data.frame(size=rep(NA, length(ls(envir = .GlobalEnv))))
  rownames(object.sizes) <- ls(envir = .GlobalEnv)
  for(i in 1:length(ls(envir = .GlobalEnv))) object.sizes$size[i] <- eval(parse(text=paste0("object.size(", ls(envir = .GlobalEnv)[i], ")")))
  object.sizes$obj <- rownames(object.sizes)
  ord <- order(object.sizes$size, decreasing = T)
  object.sizes <- object.sizes[ord,]
  return(object.sizes)
}

#' Geometric Mean
#' 
#' Function for calculating the geometric mean of a set of positive numbers.
#' @param x A vector of positive numbers.
#' @param na.rm a logical value indicating whether NA values should be stripped before the computation proceeds.
#' @return If all values of x are positive, the geometric mean is returned as a numeric vector of length one.  If 
#'    any values are non-positive, NA is returned.
#' @examples
#' x <- c(1, 10, 100)
#' mean(x)
#' geomean(x)
#' geomean(c(-1,x))
#' geomean(c(0,x))
#' geomean(c(NA,x))
#' geomean(c(NA,x), na.rm=TRUE)
#' @export
geomean <- function(x, na.rm=FALSE) {
  if (!is.numeric(x)) {
    warning("argument is not numeric: returning NA")
	return(NA)
  }
  if(na.rm) x <- x[!is.na(x)]
  if(!na.rm & sum(is.na(x)) > 0) return(NA)
  if(sum(x[!is.na(x)] <= 0) > 0) {
    warning("All values must be positive to calculate the geometric mean.")
    return(NA)
  }
  return(exp(mean(log(x))))
}

#' p-Values for \code{SSDtest} Object
#' 
#' Function for extracting the two-sided or one-sided p-values in a \code{SSDtest} object as a matrix.
#' @param x An object of class \code{SSDtest}.
#' @param alternative a character value indicating whether two-sided or one-sided p-values should be extracted. 
#'    Takes \code{"two.sided"} or \code{"one.sided"}, defaulting to \code{"two.sided"}.  If \code{"one.sided"} 
#'    then both sets of one-sided p-values are returned.
#' @return A matrix containing the corresponding p-values in \code{x}.
#' @examples
#' # Standard significance test with one fossil sample, sampling without replacement
#' test_faux_uni <- SSDtest(
#'      fossil=list("Fauxil sp. 1"=fauxil[fauxil$Species=="Fauxil sp. 1", "FHSI"]),
#'      comp=list("G. gorilla"=apelimbart[apelimbart$Species=="Gorilla gorilla", "FHSI"],
#'                "H. sapiens"=apelimbart[apelimbart$Species=="Homo sapiens", "FHSI"],
#'                "P. troglodytes"=apelimbart[apelimbart$Species=="Pan troglodytes", "FHSI"],
#'                "H. lar"=apelimbart[apelimbart$Species=="Hylobates lar", "FHSI"]),
#'      fossilsex=NULL,
#'      compsex=list("G. gorilla"=apelimbart[apelimbart$Species=="Gorilla gorilla", "Sex"],
#'                   "H. sapiens"=apelimbart[apelimbart$Species=="Homo sapiens", "Sex"],
#'                   "P. troglodytes"=apelimbart[apelimbart$Species=="Pan troglodytes", "Sex"],
#'                   "H. lar"=apelimbart[apelimbart$Species=="Hylobates lar", "Sex"]),
#'      methsUni=c("MMR", "BDI"),
#'      limit=1000,
#'      nResamp=100)
#' pvals(test_faux_uni) # defaults to two-sided
#' pvals(test_faux_uni, alternative="one.sided")
#' @export
pvals <- function(x, alternative="two.sided") {
  if (!inherits(x, "SSDtest")) stop("'x' must be an object of class 'SSDtest'.")
  alternative <- match.arg(alternative, choices=c("two.sided", "one.sided"))
  if(alternative=="two.sided") {
    out <- x$pvalues$p.twosided
	#rownames(out) <- paste0("H0: ", rownames(out), " = 0")
	rownames(out) <- gsub("-", "=", rownames(out))
	rownames(out) <- paste0("H0: ", rownames(out))
  }
  if(alternative=="one.sided") {
    out <- x$pvalues$p.onesided
	#rownames(out) <- paste0("H0: ", rownames(out), " <= 0")
	rownames(out) <- gsub("-", "<=", rownames(out))
	rownames(out) <- paste0("H0: ", rownames(out))
  }
  if (!is.null(attr(x, "warntext"))) warning(attr(x, "warntext"))
  return(out)
}

#' Central Tendencies for Resampled Distributions in \code{SSDtest} Object
#' 
#' Function for extracting the means and medians of resampled values from a \code{SSDtest} object as a matrix.
#' @param x An object of class \code{SSDtest}.
#' @return A matrix containing the means and medians of resampled distributions in \code{x}.
#' @examples
#' # Standard significance test with one fossil sample, sampling without replacement
#' test_faux_uni <- SSDtest(
#'      fossil=list("Fauxil sp. 1"=fauxil[fauxil$Species=="Fauxil sp. 1", "FHSI"]),
#'      comp=list("G. gorilla"=apelimbart[apelimbart$Species=="Gorilla gorilla", "FHSI"],
#'                "H. sapiens"=apelimbart[apelimbart$Species=="Homo sapiens", "FHSI"],
#'                "P. troglodytes"=apelimbart[apelimbart$Species=="Pan troglodytes", "FHSI"],
#'                "H. lar"=apelimbart[apelimbart$Species=="Hylobates lar", "FHSI"]),
#'      fossilsex=NULL,
#'      compsex=list("G. gorilla"=apelimbart[apelimbart$Species=="Gorilla gorilla", "Sex"],
#'                   "H. sapiens"=apelimbart[apelimbart$Species=="Homo sapiens", "Sex"],
#'                   "P. troglodytes"=apelimbart[apelimbart$Species=="Pan troglodytes", "Sex"],
#'                   "H. lar"=apelimbart[apelimbart$Species=="Hylobates lar", "Sex"]),
#'      methsUni=c("MMR", "BDI"),
#'      limit=1000,
#'      nResamp=100)
#' centers(test_faux_uni)
#' @export
centers <- function(x) {
  if (!inherits(x, "SSDtest")) stop("'x' must be an object of class 'SSDtest'.")
  return(x$methcombos)
}

#' @export
confint.dimorphResampledUni <- function(x, conf.level=0.95, alternative="two.sided", type="estimate") {
  type <- match.arg(type, choices=c("estimate", "bias"))
  alternative <- match.arg(alternative, choices=c("two.sided", "less", "greater"))
  if (class(conf.level)!="numeric") stop("'conf.level' must be a number between 0 and 1.")
  if (conf.level <= 0 | conf.level >= 1) stop("'conf.level' must be a number between 0 and 1.")
  if (!is.null(attr(x, "resampling"))) if (attr(x, "resampling")!="bootstrap") warning("This object was not generated by a bootstrap -\ntake care when interpreting these intervals.")
  if(type=="estimate") {
    #if (is.null(x$CI)) stop("There are no confidence intervals associated with this object.")
    {if (!is.null(x$CI) & attr(x, "conf.level")==conf.level & attr(x, "alternative")==alternative) {
	  out <- x$CI
	}
	else {
	  if (!is.null(x$CI) & attr(x, "conf.level")!=conf.level) warning("These intervals are for a different confidence level than\noriginally calculated when the bootstrap was run.")
	  if (!is.null(x$CI) & attr(x, "alternative")!=alternative) warning("These intervals are for a different value of 'alternative' than\noriginally calculated when the bootstrap was run.")
	  out <- tapply(x$estimates$estimate, INDEX=x$estimates$methodUni:x$estimates$center,
	                FUN=dimorph:::getCI, alternative=alternative, conf.level=conf.level)
	  out <- do.call(rbind, out)
	  out <- out[,c("lower_lim", "upper_lim")]
	  nms <- do.call(rbind, strsplit(rownames(out), split=":"))
	  out <- data.frame(methodUni=nms[,1], center=nms[,2], out)
	  out$methodUni <- factor(out$methodUni, levels=levels(x$estimates$methodUni))
	  out$center <- factor(out$center, levels=levels(x$estimates$center))
	  rownames(out) <- NULL
	}} # end if-else
  }
  if(type=="bias") {
    #if (is.null(x$CIbias)) stop("There are no bias confidence intervals associated with this object.")
    {if (!is.null(x$CIbias) & attr(x, "conf.level")==conf.level & attr(x, "alternative")==alternative) {
	  out <- x$CIbias
	}
	else {
	  if (!is.null(x$CIbias) & attr(x, "conf.level")!=conf.level) warning("These intervals are for a different confidence level than\noriginally calculated when the bootstrap was run.")
	  if (!is.null(x$CIbias) & attr(x, "alternative")!=alternative) warning("These intervals are for a different value of 'alternative' than\noriginally calculated when the bootstrap was run.")
	  out <- tapply(x$estimates$bias, INDEX=x$estimates$methodUni:x$estimates$center,
	                FUN=dimorph:::getCI, alternative=alternative, conf.level=conf.level)
	  out <- do.call(rbind, out)
	  out <- out[,c("lower_lim", "upper_lim")]
	  nms <- do.call(rbind, strsplit(rownames(out), split=":"))
	  out <- data.frame(methodUni=nms[,1], center=nms[,2], out)
	  out$methodUni <- factor(out$methodUni, levels=levels(x$estimates$methodUni))
	  out$center <- factor(out$center, levels=levels(x$estimates$center))
	  out <- out[!(out$methodUni %in% c("SSD", "CVsex", "sdlogsex")),]
	  rownames(out) <- NULL
	}} # end if-else
  }
  colnames(out) <- gsub("_lim", paste0("_lim", conf.level), colnames(out))
  return(out)
}

#' @export
confint.dimorphResampledMulti <- function(x, conf.level=0.95, alternative="two.sided", type="estimate") {
  type <- match.arg(type, choices=c("estimate", "bias"))
  alternative <- match.arg(alternative, choices=c("two.sided", "less", "greater"))
  if (class(conf.level)!="numeric") stop("'conf.level' must be a number between 0 and 1.")
  if (conf.level <= 0 | conf.level >= 1) stop("'conf.level' must be a number between 0 and 1.")
  if (!is.null(attr(x, "resampling"))) if (attr(x, "resampling")!="bootstrap") warning("This object was not generated by a bootstrap -\ntake care when interpreting these intervals.")
  if(type=="estimate") {
    #if (is.null(x$CI)) stop("There are no confidence intervals associated with this object.")
    {if (!is.null(x$CI) & attr(x, "conf.level")==conf.level & attr(x, "alternative")==alternative) {
	  out <- x$CI
	}
	else {
	  if (!is.null(x$CI) & attr(x, "conf.level")!=conf.level) warning("These intervals are for a different confidence level than\noriginally calculated when the bootstrap was run.")
	  if (!is.null(x$CI) & attr(x, "alternative")!=alternative) warning("These intervals are for a different value of 'alternative' than\noriginally calculated when the bootstrap was run.")
	  out <- tapply(x$estimates$estimate,
	                INDEX=x$estimates$methodUni:x$estimates$methodMulti:x$estimates$center:x$estimates$datastructure,
				    FUN=dimorph:::getCI, alternative=alternative, conf.level=conf.level)
	  out <- do.call(rbind, out)
	  #out <- out[,c("lower_lim", "upper_lim")]
	  nms <- do.call(rbind, strsplit(rownames(out), split=":"))
	  out <- data.frame(methodUni=nms[,1], methodMulti=nms[,2], center=nms[,3], datastructure=nms[,4],
	                    proportionNA=out[,"proportionNA"], lower_lim=out[,"lower_lim"], upper_lim=out[,"upper_lim"])
	  out$methodUni <- factor(out$methodUni, levels=levels(x$estimates$methodUni))
	  out$methodMulti <- factor(out$methodMulti, levels=levels(x$estimates$methodMulti))
	  out$center <- factor(out$center, levels=levels(x$estimates$center))
	  out$datastructure <- factor(out$datastructure, levels=levels(x$estimates$datastructure))
	  rownames(out) <- NULL
	}} # end if-else
  }
  if(type=="bias") {
    #if (is.null(x$CIbias)) stop("There are no bias confidence intervals associated with this object.")
    {if (!is.null(x$CIbias) & attr(x, "conf.level")==conf.level & attr(x, "alternative")==alternative) {
	  out <- x$CIbias
	}
	else {
	  if (!is.null(x$CIbias) & attr(x, "conf.level")!=conf.level) warning("These intervals are for a different confidence level than\noriginally calculated when the bootstrap was run.")
	  if (!is.null(x$CIbias) & attr(x, "alternative")!=alternative) warning("These intervals are for a different value of 'alternative' than\noriginally calculated when the bootstrap was run.")
	  out <- tapply(x$estimates$bias,
	                INDEX=x$estimates$methodUni:x$estimates$methodMulti:x$estimates$center:x$estimates$datastructure,
				    FUN=dimorph:::getCI, alternative=alternative, conf.level=conf.level)
	  out <- do.call(rbind, out)
	  #out <- out[,c("lower_lim", "upper_lim")]
	  nms <- do.call(rbind, strsplit(rownames(out), split=":"))
	  out <- data.frame(methodUni=nms[,1], methodMulti=nms[,2], center=nms[,3], datastructure=nms[,4],
	                    proportionNA=out[,"proportionNA"], lower_lim=out[,"lower_lim"], upper_lim=out[,"upper_lim"])
	  out$methodUni <- factor(out$methodUni, levels=levels(x$estimates$methodUni))
	  out$methodMulti <- factor(out$methodMulti, levels=levels(x$estimates$methodMulti))
	  out$center <- factor(out$center, levels=levels(x$estimates$center))
	  out$datastructure <- factor(out$datastructure, levels=levels(x$estimates$datastructure))
	  out <- out[!(out$methodUni %in% c("SSD", "CVsex", "sdlogsex")),]
	  rownames(out) <- NULL
	}} # end if-else
  }
  colnames(out) <- gsub("_lim", paste0("_lim", conf.level), colnames(out))
  return(out)
}

