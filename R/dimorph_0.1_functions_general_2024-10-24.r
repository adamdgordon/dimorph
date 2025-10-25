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

#' @export
centers <- function(x) {
  if (!inherits(x, "SSDtest")) stop("'x' must be an object of class 'SSDtest'.")
  return(x$methcombos)
}

#' @export
confint.dimorphResampledUni <- function(x, type="estimate") {
  type <- match.arg(type, choices=c("estimate", "bias"))
  if(type=="estimate") {
    if (is.null(x$CI)) stop("There are no confidence intervals associated with this object.")
	out <- x$CI
  }
  if(type=="bias") {
    if (is.null(x$CIbias)) stop("There are no bias confidence intervals associated with this object.")
    out <- x$CIbias
  }
  colnames(out) <- gsub("_lim", paste0("_lim", attr(x, "conf.level")), colnames(out))
  return(out)
}

#' @export
confint.dimorphResampledMulti <- function(x, type="estimate") {
  type <- match.arg(type, choices=c("estimate", "bias"))
  if(type=="estimate") {
    if (is.null(x$CI)) stop("There are no confidence intervals associated with this object.")
	out <- x$CI
  }
  if(type=="bias") {
    if (is.null(x$CIbias)) stop("There are no bias confidence intervals associated with this object.")
    out <- x$CIbias
  }
  colnames(out) <- gsub("_lim", paste0("_lim", attr(x, "conf.level")), colnames(out))
  return(out)
}

