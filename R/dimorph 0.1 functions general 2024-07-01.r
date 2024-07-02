#' @export
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
#' @return If all values of x are positive, the geometric mean is returned as a numeric vector of length one.  If any values are non-positive, NA is returned.
#' @examples
#' x <- c(1, 10, 100)
#' mean(x)
#' geomean(x)
#' geomean(c(-1,x))
#' geomean(c(0,x))
#' geomean(c(NA,x))
#' @export
geomean <- function(x, na.rm=F) {
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

