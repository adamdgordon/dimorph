% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dimorph_0.1_functions_analyses_2025-06-07.r
\name{resampleSSD}
\alias{resampleSSD}
\title{Resample Univariate or Multivariate Dimorphism}
\usage{
resampleSSD(
  x,
  struc = NULL,
  compsex = NULL,
  npersample = NA,
  nResamp = 10000,
  exact = F,
  limit = 5e+05,
  matchvars = F,
  replace = F,
  datastruc = NULL,
  methsMulti = NULL,
  methsUni = c("MMR", "BDI"),
  sex.female = 1,
  center = "geomean",
  templatevar = NULL,
  na.rm = T,
  ncorrection = F,
  details = F
)
}
\arguments{
\item{x}{A matrix or data frame of measurements from a comparative sample, with rows
corresponding to individual specimens and columns corresponding to size variables.  Sex data
should not be included.  Should not include \code{NA}s; resampling addresses will not be
affected, but these addresses may generate errors when used by other functions.}

\item{struc}{Structure information for the data set being compared against, typically from one or 
more fossil samples with missing data. \code{struc} must be either (1) a matrix or dataframe
of measurements (which can include \code{NA}s) or (2) a vector of integer sample sizes for each
variable.}

\item{compsex}{A vector indicating sex for the individuals in \code{x}.  Sex
information is not included in any calculations in this function but will be included
as a list element in the returned object.  Defaults to \code{NULL}.}

\item{npersample}{Integer specifying the sample size of resampled datasets.  Defaults to \code{NA}.
If set to \code{NA}, the sample size of resampled datasets is set equal to the sample size of 
\code{x}.}

\item{nResamp}{Integer specifying the number of resampling iterations to calculate addresses
for if Monte Carlo sampling is used.}

\item{exact}{Logical scalar specifying whether to sample all unique combinations of sample
size \code{npersample} from \code{x}.  Defaults to \code{FALSE}.  If set to \code{FALSE}, or if set 
to \code{TRUE} and the number of unique combinations exceeds \code{limit}, then Monte Carlo sampling
is used instead.}

\item{limit}{Integer setting the upper limit on the number of unique combinations allowable 
for exact resampling.  If exact resampling would produce more resampled datasets than this number, 
Monte Carlo resampling is used instead.  Defaults to 500,000.}

\item{matchvars}{Logical scalar specifying whether to compare the variable names in
\code{comparative} and \code{struc} and pare them both down to the set of shared variables. 
If \code{FALSE}  and variable names differ then an error will be returned.  Defaults to
\code{FALSE}.}

\item{replace}{Logical scalar passed to \code{\link[base]{sample}} specifying whether or not
to sample with replacement.  Defaults to \code{FALSE}.}

\item{datastruc}{If multivariate data are used, this is a character string specifiying whether to 
incorporate the missing data structure
into dimorphism estimates (\code{"missing"}), whether to downsample to the missing data sample size but 
keep all metric data for the comparative sample (\code{"complete"}), or to perform both types of
resampling separately (\code{"both"}).  Defaults to \code{"missing"}, and ignored if only univariate
data are provided.}

\item{methsMulti}{A character string specifying the multivariate method used to 
calculate or estimate dimorphism.  Note that regardless of the value of this argument,
multivariate estimation procedures will only be carried out if \code{x} is a multivariate 
dataset.  See \code{\link{dimorph}} for options.}

\item{methsUni}{A character string specifying the univariate method used to calculate or estimate dimorphism. 
See \code{\link{dimorph}} for options.}

\item{sex.female}{An integer scalar (1 or 2) specifying which level of \code{sex} 
corresponds to female.  Ignored if \code{sex} is \code{NULL}.  Defaults to 1.}

\item{center}{A character string specifying the method used to calculate a mean, either \code{"geomean"} 
(default) which uses the geometric mean, or \code{"mean"} which uses the arithmetic mean.  More broadly, 
\code{"geomean"} indicates analyses are conducted in logarithmic data space and \code{"mean"} indicates 
analyses are conducted in raw data space.  Some methods can only be applied in one domain or the other: 
\code{"CV"} and \code{"CVsex"} are always calculated in raw data space and \code{center} will be set to 
\code{"mean"} for these methods regardless of the value set by the user; \code{"MoM"}, \code{"sdlog"}, 
and \code{"sdlogsex"} are always calculated in logarithmic data space and \code{center} will be set to 
\code{"geomean"} for these methods regardless of the value set by the user.}

\item{templatevar}{A character object or integer value specifying the name or column number 
of the variable in \code{x} to be estimated using the template method.  Ignored if template method 
is not used.  Defaults to \code{NULL}.}

\item{na.rm}{A logical scalar indicating whether NA values should be stripped before
the computation proceeds.  Defaults to \code{TRUE}.}

\item{ncorrection}{A logical scalar indicating whether to apply Sokal and Braumann's (1980) 
size correction factor to CV estimates.  Defaults to \code{FALSE}.}

\item{details}{A logical scalar indicating whether variable name and specimen names 
should be retained (if available) as attributes in the output object.  Defaults to 
\code{FALSE}.}
}
\value{
A list of class \code{dimorphResampledUni} or \code{dimorphResampledMulti} containing a dataframe
  with resampled dimorphism estimates and a \code{dimorphAds} object containg resampled addresses produced
  by \code{\link{getsampleaddresses}}.  Plotting this object produces violin plots for
  all resampled distributions.
}
\description{
Function to generate a set of resampled dimorphism estimates for a univariate or multivariate sample.
  Called by \code{\link{SSDtest}}.
}
\examples{
## Univariate
data(apelimbart)
gor <- apelimbart[apelimbart$Species=="Gorilla gorilla",]
# this is effectively a bootstrap, although see 'bootdimorph'
gorSSD <- resampleSSD(gor[,"FHSI", drop=FALSE], methsUni=c("SSD", "MMR", "BDI"),
                      compsex=gor$Sex, nResamp=100, replace=TRUE)
gorSSD
plot(gorSSD)

# now downsample to fossil sample size and sample without replacement
SSDvars <- c("HHMaj")
gorSSD1 <- resampleSSD(x=apelimbart[apelimbart$Species=="Gorilla gorilla", SSDvars, drop=FALSE],
                       struc=fauxil[fauxil$Species=="Fauxil sp. 1", SSDvars, drop=FALSE],
                       compsex=apelimbart[apelimbart$Species=="Gorilla gorilla", "Sex"],
                       exact=TRUE, matchvars=TRUE, replace=FALSE, methsUni=c("SSD", "MMR", "BDI"))
gorSSD1
plot(gorSSD1)

# or run 'getsampleaddresses' first
addressesUni <- getsampleaddresses(
     comparative=apelimbart[apelimbart$Species=="Gorilla gorilla", SSDvars, drop=FALSE],
     struc=fauxil[fauxil$Species=="Fauxil sp. 1", SSDvars, drop=FALSE],
     compsex=apelimbart[apelimbart$Species=="Gorilla gorilla", "Sex"],
     exact=TRUE, matchvars=TRUE, replace=FALSE)
gorSSD2 <- resampleSSD(addressesUni, methsUni=c("SSD", "MMR", "BDI"))
gorSSD2
plot(gorSSD2)

## Multivariate
SSDvars <- c("HHMaj","RHMaj","FHSI","TPML")
gorSSDmulti1 <- resampleSSD(x=apelimbart[apelimbart$Species=="Gorilla gorilla", SSDvars],
                            struc=fauxil[fauxil$Species=="Fauxil sp. 1", SSDvars],
                            compsex=apelimbart[apelimbart$Species=="Gorilla gorilla", "Sex"],
                            nResamp=100,
                            datastruc="both",
                            methsMulti = c("GMM"),
                            methsUni = c("SSD", "MMR", "BDI"),
                            matchvars=TRUE,
                            replace=FALSE)
gorSSDmulti1
plot(gorSSDmulti1)

# or run 'getsampleaddresses' first
addresses <- getsampleaddresses(comparative=
               apelimbart[apelimbart$Species=="Gorilla gorilla", SSDvars],
               struc=fauxil[fauxil$Species=="Fauxil sp. 1", SSDvars],
               compsex=apelimbart[apelimbart$Species=="Gorilla gorilla", "Sex"],
               nResamp=100, matchvars=TRUE, replace=FALSE)
gorSSDmulti2 <- resampleSSD(x=addresses,
                            datastruc="both",
                            methsMulti = c("GMM"),
                            methsUni = c("SSD", "MMR", "BDI"))
gorSSDmulti2
plot(gorSSDmulti2)
}
\seealso{
\code{\link{bootdimorph}}, \code{\link{dimorph}}, \code{\link{getsampleaddresses}}
}
