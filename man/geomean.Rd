% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dimorph_0.1_functions_general_2024-07-16.r
\name{geomean}
\alias{geomean}
\title{Geometric Mean}
\usage{
geomean(x, na.rm = FALSE)
}
\arguments{
\item{x}{A vector of positive numbers.}

\item{na.rm}{a logical value indicating whether NA values should be stripped before the computation proceeds.}
}
\value{
If all values of x are positive, the geometric mean is returned as a numeric vector of length one.  If 
   any values are non-positive, NA is returned.
}
\description{
Function for calculating the geometric mean of a set of positive numbers.
}
\examples{
x <- c(1, 10, 100)
mean(x)
geomean(x)
geomean(c(-1,x))
geomean(c(0,x))
geomean(c(NA,x))
geomean(c(NA,x), na.rm=TRUE)
}
