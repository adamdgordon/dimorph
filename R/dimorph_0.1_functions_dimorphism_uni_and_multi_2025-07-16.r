#' Estimate Dimorphism in a Univariate or Multivariate Sample
#' 
#' Function to calculate or estimate dimorphism for a univariate or multivariate sample.
#' @param x A dataframe, matrix, or vector of positive numbers corresponding to measurements for one or more size variable(s).  Values must be in the original measurement space, \bold{not} log-transformed.
#' @param method A character string specifying the univariate method used to calculate or estimate dimorphism.  Options 
#'   include:
#'   \itemize{
#'   \item \code{"SSD"}: \bold{Sexual Size Dimorphism}. Follows Smith (1999). Calculates actual sexual 
#'     dimorphism in the sample as the ratio of
#'     mean male size to mean female size.  Depending on \code{center}, sex-specific means are calculated
#'     either as geometric means or arithmetic means.  Requires \code{sex} to be specified.
#'   \item \code{"MMR"}: \bold{Mean Method Ratio}. Follows Godfrey et al. (1993). Splits the sample at its mean, 
#'     then calculates the ratio of the mean 
#'     of measurements larger than the overall mean to the mean of measurements smaller than the
#'     overall mean.  If any measurements are exactly equal to the overall mean, they contribute to both
#'     the larger and smaller group as half an individual in a weighted mean.  Depending on \code{center}, 
#'     the overall mean and subgroup means are calculated either as geometric means or arithmetic 
#'     means.  Ignores \code{sex}.
#'   \item \code{"BDI"}: \bold{Binomial Dimorphism Index}. Follows Reno et al. (2003). Given \emph{n} measurements,
#'     calculates all possible ratios of the mean of larger specimens to the mean of
#'     smaller specimens when the sample is split into the \emph{k} largest specimens and
#'     \emph{n-k} smallest specimens, where \emph{k} ranges from 1 to \emph{n}-1.  A weighted 
#'     mean is then calculated for all ratios, where the weights are equal to the probability
#'     of \emph{k} successes in \emph{n} trials in a binomial distribution. Depending on 
#'     \code{center}, means are calculated either as geometric means or arithmetic 
#'     means.  Ignores \code{sex}.
#'   \item \code{"ERM"}: \bold{Exact Resampling Method}.  Modification of Lee's (2001) Assigned Resampling 
#'     Method (ARM) following Gordon (2025).  ARM is a resampling-based estimate of dimorphism that repeatedly 
#'     samples two values
#'     with replacement from \code{x}, then calculates their ratio as long as both are neither more than 0.5 standard
#'     deviations above the mean or both 0.5 standard deviations below the mean (otherwise the pair is 
#'     rejected).  ARM typically oversamples the possible combination of two values sampled from a small sample 
#'     (as originally described it samples 1,000 pairs, whereas a sample of 20 measurements only has 210 possible pairs) 
#'     and sampling with replacement biases dimorphism estimates downwards by the incorporation of multiple ratios
#'     of 1 whever the same value is sampled twice and is not rejected by retention criterion.
#'     \code{"ERM"} performs an exact resampling of all possible
#'     pairs of values without replacement, but otherwise follows Lee's algorithm.
#'     Depending on \code{center}, the procedure is applied in either logarithmic (\code{"geomean"}) or
#'     raw (\code{"mean"}) data space.  Ignores \code{sex}.
#'   \item \code{"MoM"}: \bold{Method of Moments}.  Follows Josephson et al. (1996).  Assumes that the sample is 
#'     a mixture of two underlying lognormal distributions and uses three moments around the mean of the logged
#'     combined sex distribution to estimate the means of the underlying distributions.  This calculation is 
#'     always performed on the log-transformed data regardless of the value of \code{center}.
#'     Assumes that the sample contains an equal proportion of females and males and that those two subsamples 
#'     have equal variance.  Ignores \code{sex}.
#'   \item \code{"FMA"}: \bold{Finite Mixture Analysis}. Follows Godfrey et al. (1993). Assumes that the sample is 
#'     a mixture of two underlying normal or lognormal distributions.  Assumes that the sample contains an equal
#'     proportion of females and males and that those two subsamples have equal variance, then estimates the
#'     maximum separation of the two means. Depending on \code{center}, the underlying distributions are
#'     treated as either normal (\code{"mean"}) or lognormal (\code{"geomean"}).  Ignores \code{sex}.
#'   \item \code{"BFM"}: \bold{Bayesian Finite Mixture}.  Follows Gordon (2025).  Assumes that the sample is a 
#'     finite mixture of two 
#'     underlying normal or lognormal distributions.  Unlike \code{"FMA"} and \code{"MoM"}, estimates 
#'     the proportion of females and males assuming that they may not be equal, and uses a Bayesian Information 
#'     Criterion (BIC) approach to select between a model that estimates a single variance for both sexes and a 
#'     model that estimates variances separately for the two constituent distributions using 
#'     \code{\link[mclust]{mclustBIC}}.  It then calculates the ratio of
#'     the two estimated means. Depending on \code{center}, the underlying distributions are
#'     treated as either normal (\code{"mean"}) or lognormal (\code{"geomean"}).  Ignores \code{sex}.
#'     When performed on lognormal data, this method is similar to the pdPeak method of Sasaki et al. (2021), 
#'     particularly when the BIC procedure selects an equal variance model (which it typically does).
#'   \item \code{"CV"}: \bold{Coefficient of Variation}.  Calculates the coefficient of variation as the 
#'     standard deviation of \code{x} divided by the mean of \code{x} then multiplied by 100.  This calculation
#'     is always performed on the raw data regardless of the value of \code{center} (an analogous 
#'     method using logarithmic data is \code{"sdlog"}).  Ignores \code{sex}.  Additionally, Sokal and Braumann's (1980) 
#'     size correction factor can be applied by setting \code{ncorrection} to \code{TRUE}, although this is
#'     \code{FALSE} by default.
#'   \item \code{"CVsex"}: \bold{Modified Coefficient of Variation}.  Calculates a modified 
#'     version of the coefficient of variation: the standard deviation is replaced by the 
#'     square root of the sum of squared differences of every value of \code{x} from the unweighted
#'     mean of the sex-specific means in \code{x} divided by the square root of \emph{n}-1, and this is 
#'     divided by the unweighted mean of the sex-specific means, then multiplied by 100.  This calculation
#'     is always performed on the raw data regardless of the value of \code{center} (an analogous 
#'     method using logarithmic data is \code{"sdlogsex"}).  Requires \code{sex}.  Additionally, Sokal and Braumann's (1980) 
#'     size correction factor can be applied by setting \code{ncorrection} to \code{TRUE}, although this is
#'     \code{FALSE} by default.
#'   \item \code{"sdlog"}: \bold{Standard Deviation of Logged Data}.  First, \code{x} is log-transformed
#'     using the natural logarithm, then the standard deviation is calculated.  This is a measure of
#'     proportional variation of the values of \code{x} about their geometric mean; \emph{i.e.}, analagous
#'     to the coefficient of variation for a lognormal distribution. This calculation
#'     is always performed on the log-transformed data regardless of the value of \code{center} 
#'     (an analogous method using raw data is \code{"CV"}).  Ignores \code{sex}.
#'   \item \code{"sdlogsex"} \bold{Modified Standard Deviation of Logged Data}.  First, \code{x} is 
#'     log-transformed using the natural logarithm.  Then a modified version of standard deviation is 
#'     calculated: the square root of the sum of squared differences of every logged value of \code{x} 
#'     from the unweighted mean of the sex-specific means of log-transformed \code{x}, divided by the 
#'     square root of \emph{n}-1. This calculation is always performed on the log-transformed data 
#'     regardless of the value of \code{center}  (an analogous method using raw data is \code{"CVsex"}).
#'     Requires \code{sex}.
#' } Defaults to \code{"SSD"}.
#' @param methodMulti A character string specifying the multivariate method used to 
#'   calculate or estimate dimorphism.  Note that regardless of the value of this argument,
#'   multivariate estimation procedures will only be carried out if \code{x} is a multivariate 
#'   dataset.  Options include:
#'   \itemize{
#'   \item \code{"GMM"}: Follows the Geometric Mean Method of Gordon et al. (2008).  The
#'     selected univariate method is applied to all variables, then the geometric mean
#'     is calculated of the dimorphism estimates for all variables to produce a single
#'     estimate for the whole data set. Note 
#'     that this methodology is not appropriate for variance-based univariate methods;
#'     \emph{i.e.}, \code{"CV"}, \code{"CVsex"}, \code{"sdlog"}, and \code{"sdlogsex"}.
#'   \item \code{"GMsize"}:  If \code{x} is a dataframe or matrix, this method first 
#'     calculates overall size as the geometric mean of measurements in all variables for
#'     those specimens that are complete for all variables in the data set.  The selected
#'     univariate method is then applied to this measure of overall size.
#'   \item \code{"TM"}:  Follows the template method of Reno et al. (2003). A variable of interest
#'     is specified by the user with the argument \code{templatevar}.  The algorithm identifies
#'     a template individual that can be used to estimate the largest number of values for 
#'     the selected variable of interest by using ratios between the value of that variable and
#'     other variables in the template individual, which are then multiplied by the value of those
#'     other variables in individuals missing the target variable to maximixe the dataset for 
#'     that variable.  A user-selected univariate method is then applied to the combined dataset
#'     of actual and estimated values for the target variable.
#' } Defaults to \code{"GMM"}.
#' @param sex A vector indicating sex for the measurements in \code{x}.  If present, must include 
#'   exactly two groups and have the same length as the number of specimens in \code{x}.  Non-factor 
#'   vectors will be coerced to factors if possible.  May be \code{NULL} since some methods do not 
#'   require sex information.  Methods which require sex information will generate an error if 
#'   \code{sex} is \code{NULL}.  For methods that do not require sex information, if sex is provided 
#'   it will be ignored for the calculation of the estimate, but it will be used to report the actual 
#'   proportion of females and males in the sample.  Defaults to \code{NULL}.
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
#' @param ads A vector of integer addresses for positions in the data vector \code{x} to be included 
#'   in the calculation of dimorphism; any other data in \code{x} will be ignored.  If \code{ads} is
#'   \code{NULL} then all data are included in the calculation.  Defaults to \code{NULL}.
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
#' @param dfout A logical scalar indicating whether the result should be given as a
#'   \code{dimorphEstDF} object; if \code{FALSE}, returns a \code{dimorphEst} object.  Defaults to \code{FALSE}.
#' @return Either a class \code{dimorphEst} or \code{dimorphEstDF} object. \code{dimorphEst} objects are numeric 
#'   vectors of length one corresponding to measured or estimated dimorphism in \code{x} with associated information 
#'   preserved as attributes.  \code{dimorphEstDF} objects are single-row data frames that contain the dimorphism 
#'   estimate for \code{x} along with other associated information.  Applying \code{summary} to either of these 
#'   objects provides information about the dataset and method used to generate it.
#' @seealso \code{\link[mclust]{mclustBIC}},  \code{\link[dimorph]{SSDtest}}
#' @examples
#' ## Univariate estimates:
#' data(apelimbart)
#' gorillas <- apelimbart[apelimbart$Species=="Gorilla gorilla",]
#' # Next line would generate an error: sex is required
#' # dimorph(x=gorillas$FHSI, method="SSD") 
#' gorSSD  <- dimorph(x=gorillas$FHSI, # variable and specimen names not preserved
#'                    method="SSD", sex=gorillas$Sex, details=TRUE)
#' gorSSD2 <- dimorph(x=gorillas[,"FHSI",drop=FALSE], # variable and specimen names preserved
#'                    method="SSD", sex=gorillas$Sex, details=TRUE)
#' gorSSD
#' str(gorSSD)
#' summary(gorSSD)
#' summary(gorSSD2) # results are identical to 'gorSSD'
#' summary(gorSSD, verbose=TRUE)
#' summary(gorSSD2, verbose=TRUE) # variable and specimen names preserved
#' # A subset of specimens can be specified for analysis using 'ads'
#' summary(dimorph(x=gorillas$FHSI, method="SSD", sex=gorillas$Sex, ads=c(1:10, 51:60)))
#' # Methods for estimating dimorphism:
#' summary(dimorph(x=gorillas$FHSI, method="MMR"))
#' summary(dimorph(x=gorillas$FHSI, method="BDI"))
#' summary(dimorph(x=gorillas$FHSI, method="MoM"))
#' summary(dimorph(x=gorillas$FHSI, method="FMA"))
#' summary(dimorph(x=gorillas$FHSI, method="BFM"))
#' summary(dimorph(x=gorillas$FHSI, method="ERM"))
#' summary(dimorph(x=gorillas$FHSI, method="CV"))
#' summary(dimorph(x=gorillas$FHSI, method="CV", ncorrection=TRUE))
#' summary(dimorph(x=gorillas$FHSI, method="CVsex", sex=gorillas$Sex))
#' summary(dimorph(x=gorillas$FHSI, method="CVsex", sex=gorillas$Sex, ncorrection=TRUE))
#' summary(dimorph(x=gorillas$FHSI, method="sdlog"))
#' summary(dimorph(x=gorillas$FHSI, method="sdlogsex", sex=gorillas$Sex))
#'
#' # Now setting 'dfout' to TRUE:
#' allmethods <- rbind(dimorph(x=gorillas$FHSI, method="SSD", sex=gorillas$Sex, dfout=TRUE),
#'                     dimorph(x=gorillas$FHSI, method="MMR", dfout=TRUE),
#'                     dimorph(x=gorillas$FHSI, method="BDI", dfout=TRUE),
#'                     dimorph(x=gorillas$FHSI, method="MoM", dfout=TRUE),
#'                     dimorph(x=gorillas$FHSI, method="FMA", dfout=TRUE),
#'                     dimorph(x=gorillas$FHSI, method="BFM", dfout=TRUE),
#'                     dimorph(x=gorillas$FHSI, method="ERM", dfout=TRUE),
#'                     dimorph(x=gorillas$FHSI, method="CV", dfout=TRUE),
#'                     dimorph(x=gorillas$FHSI, method="CVsex", sex=gorillas$Sex, dfout=TRUE),
#'                     dimorph(x=gorillas$FHSI, method="sdlog", dfout=TRUE),
#'                     dimorph(x=gorillas$FHSI, method="sdlogsex", sex=gorillas$Sex, dfout=TRUE))				
#' allmethods
#'
#' # Alternatively, using apply
#' res <- apply(data.frame(method=c("SSD","MMR","BDI","MoM","FMA","BFM",
#'                                  "ERM","CV","CVsex","sdlog","sdlogsex")),
#'              MARGIN=1,
#'              FUN=function(method, ...) dimorph(x=gorillas$FHSI, method=method,
#'                                                sex=gorillas$Sex, dfout=TRUE),
#'              simplify=FALSE)
#' as.data.frame(do.call(rbind, res))
#' 
#' ## Multivariate estimates:
#' # GMsize (only usable for complete datasets)
#' Gg.GMsize <- dimorph(x=gorillas[,c("FHSI","HHMaj","TPMAP","RHMaj")],
#'                      method="SSD", methodMulti="GMsize", sex=gorillas$Sex, details=TRUE)
#' Gg.GMsize
#' summary(Gg.GMsize)
#' 
#' # GMM (produces the same values for ratio estimators as GMsize when data are complete)
#' Gg.GMM1 <- dimorph(x=gorillas[,c("FHSI","HHMaj","TPMAP","RHMaj")],
#'                    method="SSD", methodMulti="GMM", sex=gorillas$Sex, details=TRUE)
#' # now with subset of gorilla data
#' Gg.GMM2 <- dimorph(x=gorillas[,c("FHSI","HHMaj","TPMAP","RHMaj")],
#'                    method="SSD", methodMulti="GMM", sex=gorillas$Sex,
#'                    ads=c(1:10, 51:60), details=TRUE)
#' summary(Gg.GMM1)
#' summary(Gg.GMM2)
#' 
#' ## Now with some simulated fossil data
#' SSDvars <- c("FHSI", "TPML", "TPMAP", "TPLAP", "HHMaj",
#'              "HHMin", "RHMaj", "RHMin", "RDAP", "RDML")
#' Fs1 <- fauxil[fauxil$Species=="Fauxil sp. 1", SSDvars]
#' Fs1GMM <- dimorph(x=Fs1, method="MMR", methodMulti="GMM", details=TRUE)
#' Fs1TM <- dimorph(x=Fs1, method="MMR", methodMulti="TM", 
#'                  templatevar="FHSI", details=TRUE)
#' Fs1GMM
#' Fs1TM
#' summary(Fs1GMM)
#' summary(Fs1TM, verbose=TRUE)
#' 
#' Fs1Both <- rbind(dimorph(x=Fs1, method="MMR", methodMulti="GMM", details=TRUE, dfout=TRUE),
#'                  dimorph(x=Fs1, method="MMR", methodMulti="TM", 
#'                          templatevar="FHSI", details=TRUE, dfout=TRUE))
#' Fs1Both
#' @references Godfrey LR, Lyon SK, Sutherland MR. (1993) Sexual dimorphism in large-bodied 
#'   primates: The case of the subfossil lemurs. \emph{American Journal of Physical 
#'   Anthropology}. 90:315-334. (\href{https://doi.org/10.1002/ajpa.1330900306}{https://doi.org/10.1002/ajpa.1330900306})
#' @references Gordon AD. (2025a) Interpreting statistical significance in hominin dimorphism: Power and Type I error 
#'   rates for resampling tests of univariate and missing-data multivariate size dimorphism estimation methods in the 
#'   fossil record. \emph{Journal of Human Evolution}. 199:103630. 
#'   (\href{https://doi.org/10.1016/j.jhevol.2024.103630}{https://doi.org/10.1016/j.jhevol.2024.103630})
#' @references Gordon AD, Green DJ, Richmond BG. (2008) Strong postcranial size dimorphism in 
#'   \emph{Australopithecus afarensis}: Results from two new resampling methods for multivariate 
#'   data sets with missing data. \emph{American Journal of Physical Anthropology}. 135:311-328.
#'   (\href{https://doi.org/10.1002/ajpa.20745}{https://doi.org/10.1002/ajpa.20745}) 
#' @references Josephson SC, Juell KE, Rogers AR. (1996) Estimating sexual dimorphism by 
#'   method-of-moments. \emph{American Journal of Physical Anthropology}. 100:191-206. 
#'   (\href{https://doi.org/10.1002/(SICI)1096-8644(199606)100:2<191::AID-AJPA3>3.0.CO;2-0}{https://doi.org/10.1002/(SICI)1096-8644(199606)100:2<191::AID-AJPA3>3.0.CO;2-0})
#' @references Lee S-H. (2001) Assigned Resampling Method: A new method to estimate size sexual dimorphism in 
#'   samples of unknown sex. \emph{Anthropological Review}. 64:21–39. 
#'   (\href{https://doi.org/10.18778/1898-6773.64.02}{https://doi.org/10.18778/1898-6773.64.02})
#' @references Reno PL, Meindl RS, McCollum MA, Lovejoy CO. (2003) Sexual dimorphism in 
#'   \emph{Australopithecus afarensis} was similar to that of modern humans. \emph{Proceedings 
#'   of the National Academy of Sciences}. 100:9404-9409. 
#'   (\href{https://doi.org/10.1073/pnas.1133180100}{https://doi.org/10.1073/pnas.1133180100})
#' @references Sasaki T, Semaw S, Rogers MJ, Simpson SW, Beyene Y, Asfaw B, White TD, Suwa G. (2021) Estimating sexual 
#'   size dimorphism in fossil species from posterior probability densities. \emph{Proceedings 
#'   of the National Academy of Sciences}. 118:e2113943118.
#'   (\href{https://doi.org/10.1073/pnas.2113943118}{https://doi.org/10.1073/pnas.2113943118})
#' @references Sokal RR, Braumann CA (1980) Significance tests for coefficients of variation and variability profiles. 
#'   \emph{Systematic Zoology}. 29:50–66. 
#'   (\href{https://doi.org/10.2307/2412626}{https://doi.org/10.2307/2412626})
#' @references Smith RJ. (1999) Statistics of sexual size dimorphism. \emph{Journal of 
#'   Human Evolution}. 36:423-458. (\href{https://doi.org/10.1006/jhev.1998.0281}{https://doi.org/10.1006/jhev.1998.0281})
#' @export
dimorph <- function(x, method="SSD", methodMulti="GMM", sex=NULL, sex.female=1,
                    center="geomean", ads=NULL, templatevar=NULL, na.rm=T, ncorrection=F,
                    details=F, dfout=F) {
  uni <- NA
  {if (class(x)[1]=="data.frame" | class(x)[1]=="matrix") {
    {if (ncol(x)==1) uni <- T
    else if (ncol(x)>1) uni <- F}
  }
  else if (class(x)[1]=="numeric" | class(x)[1]=="integer") uni <- T
  else if (class(x)[1]=="list") uni <- F}
  if (is.na(uni)) stop("'x' must be a numeric vector, matrix, dataframe, or list.")
  {if (uni) return(dimorph:::dimorphUni(x=x, method=method, sex=sex, sex.female=sex.female, center=center, ads=ads, na.rm=na.rm, ncorrection=ncorrection, details=details, dfout=dfout))
  else if (!uni) return(dimorph:::dimorphMulti(x=x, methodUni=method, methodMulti=methodMulti, sex=sex, sex.female=sex.female, center=center, ads=ads, templatevar=templatevar, na.rm=na.rm, ncorrection=ncorrection, details=details, dfout=dfout))}
  return(NA)
}


#' @noRd
dimorphUni <- function(x, method="SSD", sex=NULL, sex.female=1, center="geomean",
                       ads=NULL, na.rm=T, ncorrection=F, details=F, dfout=F) {
  method <- match.arg(method, choices=c("SSD", "MMR", "BDI", "MoM", "FMA", "BFM",
                                        "ERM", "CV", "CVsex", "sdlog", "sdlogsex"))
  center <- match.arg(center, choices=c("geomean", "mean"))
  {if (dfout) {
    out <- data.frame(estimate=NA,
	                  methodUni=method,
	                  methodMulti=NA,
	                  center=center,
	                  n.vars.overall=NA,
					  n.specimens.overall=NA,
					  proportion.female.overall=NA,
					  n.vars.realized=NA,
					  n.specimens.realized=NA,
					  proportion.female.realized=NA,
	                  proportion.missingdata.overall=NA,
					  proportion.missingdata.realized=NA,
					  proportion.templated=NA)
	attr(out$estimate, "details") <- list(ratio.means=NA,
	                                      vars.used=NA,
	                                      specimens.used=NA,
										  model.parameters=NA,
										  estvalues="raw")
    class(out) <- c("dimorphEstDF", "data.frame")
  }
  else {
    out <- NA
	attr(out, "details") <- list(estimate=NA,
	                             methodUni=method,
								 methodMulti=NA,
								 center=center,
								 n.vars.overall=NA,
								 n.specimens.overall=NA,
								 proportion.female.overall=NA,
								 n.vars.realized=NA,
								 n.specimens.realized=NA,
								 proportion.female.realized=NA,
								 proportion.missingdata.overall=NA,
								 proportion.missingdata.realized=NA,
								 proportion.templated=NA,
								 ratio.means=NA,
								 vars.used=NA,
								 specimens.used=NA,
								 model.parameters=NA,
								 estvalues="raw")
    class(out) <- c("dimorphEst", "numeric")
  }}
  varname <- NA
  specnames <- NA
  if (class(x)[1]=="data.frame") x <- as.matrix(x)
  {if (class(x)[1]=="matrix") {
    if (!inherits(x[,1], c("numeric", "integer"))) stop("All data in 'x' must be real numbers.")
    #if (!(class(x[,1])=="numeric" | class(x[,1])=="integer")) stop("All data in 'x' must be real numbers.")
	if (ncol(x) > 1) stop("dimorphUni can only accept one variable.")
	if (!is.null(colnames(x))) varname <- colnames(x)[1]
	if (is.null(rownames(x))) rownames(x) <- paste0("specimen_", 1:nrow(x))
	specnames <- rownames(x)
  }
  else if (class(x)[1]=="numeric" | class(x)[1]=="integer") {
    if (is.null(names(x))) names(x) <- paste0("specimen_", 1:length(x))
    specnames <- names(x)
  }
  else stop("'x' must be a vector or single variable dataframe with numeric or integer data.")
  }  
  if (class(x)[1]=="matrix") x <- as.vector(x)
  if(!is.null(ads)) {
    x <- x[ads]
	if(!is.null(sex)) sex <- sex[ads]
	if(!is.null(specnames)) specnames <- specnames[ads]
  }
  if (!na.rm) if (sum(is.na(x))>0) return(out)
  if (length(sort(x)) < 2) {
    warning("At least 2 values must be present.")
    return(out)
  }
  if (length(sort(x)) < 3 & method=="BFM") {
    warning("At least 3 values must be present to estimate means using 'BFM'.")
    return(out)
  }
  nspecoverall <- length(x)
  nspecrealized <- length(x)
  propFoverall <- NA
  propFrealized <- NA
  mnm <- NA
  mnf <- NA
  modpar <- NA
  {
    if(!is.null(sex)) {
      if (!inherits(sex, "factor")) sex <- droplevels(as.factor(sex))
      #if(!class(sex)=="factor") sex <- droplevels(as.factor(sex))
      if (!(nlevels(sex)==1|2)) stop("The sex variable can have no more than two levels.")
      if (nlevels(sex)==1 & method %in% c("SSD", "CVsex", "sdlogsex")) {
        warning(paste0(method, " can not be calculated when only one sex is present."))
        return(out)
      }
      if (!length(sex)==length(x)) stop("The measurement and sex vectors must have equal length.")
      if (!(sex.female==1|sex.female==2)) stop("sex.female must be either 1 or 2.")
	  propFoverall <- sum(sex==levels(sex)[sex.female]) / nspecoverall
      if (na.rm) {
        y <- data.frame(x=x, sex=sex)
		cmplt <- stats::complete.cases(y)
        y <- y[cmplt,]
        if (nrow(y) < 2) return(NA)
        x <- y$x
        sex <- y$sex
		if(!is.null(specnames)) specnames <- specnames[cmplt]
        rm(y)
		nspecrealized <- length(x)
      }
	  propFrealized <- sum(sex==levels(sex)[sex.female])/length(sex)
    }
    else {
	  cmplt <- !is.na(x)
	  x <- x[cmplt]
	  if(!is.null(specnames)) specnames <- specnames[cmplt]
	  nspecrealized <- length(x)
	}
  }
  if(method %in% c("MMR", "BDI", "MoM", "FMA", "BFM", "ERM", "CV", "sdlog")) {
    sex <- NULL
  }
  {if(center=="geomean") {fcenter=geomean}
    else if (center=="mean") {fcenter=mean}
  }
  {# start if loops for specific estimator function
    if (method=="SSD") { # start SSD
      if (is.null(sex)) stop("A variable identifying sex in specimens must be included.")
	  mnm <- fcenter(x[sex==levels(sex)[-sex.female]])
	  mnf <- fcenter(x[sex==levels(sex)[sex.female]])
      est <- mnm / mnf
    } # end SSD
    else if (method=="MMR") { # start MMR
      m <- round(fcenter(x, na.rm=na.rm),10)
      ncent <- sum(x==m, na.rm=na.rm)
      {if (ncent > 0) { # any values exactly at the mean count half towards both groups
        m1 <- x[x>m]
        m2 <- x[x<m]
        if(center=="geomean") {
          m1 <- log(m1)
          m2 <- log(m2)
          m <- log(m)
		  mnm <- exp((sum(m1)+ncent/2*m)/(length(m1)+ncent/2))
		  mnf <- exp((sum(m2)+ncent/2*m)/(length(m2)+ncent/2))
          est <- mnm / mnf
        }
        else {
		  mnm <- ((sum(m1)+ncent/2*m)/(length(m1)+ncent/2))
		  mnf <- ((sum(m2)+ncent/2*m)/(length(m2)+ncent/2))
          est <- mnm / mnf
        }
      }
        else {
		  mnm <- fcenter(x[x>=m])
		  mnf <- fcenter(x[x<m])
		  est <- mnm / mnf
		}
      }
    } # end MMR
    else if (method=="BDI") { # start BDI
      n <- length(x)
      w <- dbinom(1:(n-1), n, 0.5)
      r <- matrix(0, n-1)
      y <- sort(x)
      {
        if (center=="geomean") {
          y <- log(y)
          for (i in 1:(n-1)) {
            r[i] <- mean(y[(i+1):n]) - mean(y[1:i])
          }
          est <- exp(sum(r*w)/sum(w))
        }
        else if (center=="mean") {
          for (i in 1:(n-1)) {
            r[i] <- fcenter(y[(i+1):n]) / fcenter(y[1:i])
          }
          est <- (sum(r*w)/sum(w))
        }
      }
    } # end BDI
    else if (method=="ERM") { # start ERM - a modification of Lee's ARM
      xp <- x
      if (center=="geomean") xp <- log(xp)
      lims <- c(mean(xp)-sd(xp)*0.5, mean(xp)+sd(xp)*0.5)
      rjct <- function(z) prod(z < lims[1]) | prod(z > lims[2])
      y <- combn(xp, 2)
      y <- y[,!apply(y, 2, rjct)]
      {
        if (center=="geomean") {div <- function(z) max(z)-min(z)}
        else if (center=="mean") {div <- function(z) max(z)/min(z)}
      }
      if (length(y)==0) return(est)
      {
        if (length(y)==2) est <- div(y)
        else est <- mean(apply(y, 2, div))
      }
      if (center=="geomean") est <- exp(est)
    } # end ERM
    else if (method=="MoM") { # start MoM
      center <- "geomean"
      fcenter <- geomean
      y  <- log(x)
      z  <- (y-mean(y))/sd(y)
      m4 <- mean(z^4)
      D  <- 1.5 - m4/2
      {
        if (is.na(D)) est <- NA
        else {
          if (D < 0) D <- 0
          est <- exp(2*sd(y)*D^0.25)
        }
      }
    } # end MoM
    else if (method=="FMA") { # start FMA; Godfrey et al., 1993
      y <- x
      if (center=="geomean") y <- log(y)
      n <- length(y)
      {if (n > 100) est <- NA
	  else {
        R <- max(y) - min(y)
        {
          if (R==0) est <- NA
          else {
            k <- dimorph::FMAtable$k[dimorph::FMAtable$n==n]
            #k <- FMAtable$k[FMAtable$n==n]
            sigma <- R/(k*sqrt(2))
            {
              if (center=="geomean") {
			    mnm <- exp(mean(y)+sigma)
			    mnf <- exp(mean(y)-sigma)
			    est <- mnm / mnf
			  }
              else if (center=="mean") {
			    mnm <- mean(y)+sigma
			    mnf <- mean(y)-sigma
			    est <- mnm / mnf
			  }
            }
          }
        }
	  }}
    } # end FMA
    else if (method=="BFM") { # start BFM
      y <- x
      if (center=="geomean") y <- log(y)
      VM <- mclust::mclustBIC(data = y, ModelNames="V", verbose=F)
      Gaussian.two <- mclust::summaryMclustBIC(VM, data = y, G=2)
	  {if (is.null(Gaussian.two$parameters)) est <- NA
	  else {
      {
        if (center=="geomean") {
		  mnm <- exp(max(Gaussian.two$parameters$mean))
		  mnf <- exp(min(Gaussian.two$parameters$mean))
		  est <- mnm / mnf
		}
        else if (center=="mean") {
		  mnm <- max(Gaussian.two$parameters$mean)
		  mnf <- min(Gaussian.two$parameters$mean)
		  est <- mnm / mnf
		}
      }
	  modpar <- Gaussian.two$parameters
	  if (modpar$mean[1] > modpar$mean[2]) {
	    modpar$mean <- rev(modpar$mean)
	    modpar$pro <- rev(modpar$pro)
		if (modpar$variance$modelName=="V") {
		  modpar$variance$sigmasq <- rev(modpar$variance$sigmasq)
		  modpar$variance$scale <- rev(modpar$variance$scale)
		}
	  }}}
    } # end BFM
    else if (method=="CV") { # start CV
      center="mean"
      fcenter=mean
      est <- sd(x)/mean(x)*100
	  if (ncorrection) est <- est*(1+1/(4*length(x)))
    } # end CV
    else if (method=="CVsex") { # start CVsex
      center="mean"
      fcenter=mean
      mid <- mean(c(mean(x[sex==levels(sex)[1]]), mean(x[sex==levels(sex)[2]])))
      est <- sqrt(sum((x-mid)^2)/(length(x)-1))/mid*100
 	  if (ncorrection) est <- est*(1+1/(4*length(x)))
   } # end CVsex
    else if (method=="sdlog") { # start sdlog
      center="geomean"
      fcenter=geomean
      y <- log(x)
      est <- sd(y)
    } # end sdlog
    else if (method=="sdlogsex") { # start sdlogsex
      center="geomean"
      fcenter=geomean
      y <- log(x)
      mid <- mean(c(mean(y[sex==levels(sex)[1]]), mean(y[sex==levels(sex)[2]])))
      est <- sqrt(sum((y-mid)^2)/(length(y)-1))
    } # end sdlogsex
  }# end if loops for specific estimator function
  {if (!dfout) {
    out[[1]] <- est
	names(out) <- method
    attr(out, "details")$methodUni <- method
    attr(out, "details")$methodMulti <- NA
    attr(out, "details")$center <- center
    attr(out, "details")$n.vars.overall <- 1
    attr(out, "details")$n.specimens.overall <- nspecoverall
    attr(out, "details")$proportion.female.overall <- propFoverall
    attr(out, "details")$n.vars.realized <- 1
    attr(out, "details")$n.specimens.realized <- length(x)
    attr(out, "details")$proportion.female.realized <- propFrealized
    attr(out, "details")$proportion.missingdata.overall <- 1-(length(x)/nspecoverall)
    attr(out, "details")$proportion.missingdata.realized <- 0
    attr(out, "details")$proportion.templated <- NA
	attr(out, "details")$ratio.means <- c(numerator=mnm, denominator=mnf)
	attr(out, "details")$vars.used <- NA
	attr(out, "details")$specimens.used <- NA
	attr(out, "details")$model.parameters <- modpar	
	if (method=="CV" | method=="CVsex") attr(out, "details")$ncorrection <- ncorrection	
    if(details) {
      attr(out, "details")$vars.used <- varname
	  attr(out, "details")$specimens.used <- specnames
    }      
    #class(out) <- c("dimorphEst", "numeric")
    #class(out) <- "dimorphEst"
    }
  else {
    out$estimate[[1]] <- est
	rownames(out) <- method
    out$methodUni <- method
    out$methodMulti <- NA
    out$center <- center
    out$n.vars.overall <- 1
    out$n.specimens.overall <- nspecoverall
    out$proportion.female.overall <- propFoverall
    out$n.vars.realized <- 1
    out$n.specimens.realized <- length(x)
    out$proportion.female.realized <- propFrealized
    out$proportion.missingdata.overall <- 1-(length(x)/nspecoverall)
    out$proportion.missingdata.realized <- 0
    out$proportion.templated <- NA
	attr(out$estimate, "details")$ratio.means <- c(numerator=mnm, denominator=mnf)
	attr(out$estimate, "details")$vars.used <- NA
	attr(out$estimate, "details")$specimens.used <- NA
	attr(out$estimate, "details")$model.parameters <- modpar	
	if (method=="CV" | method=="CVsex") attr(out, "details")$ncorrection <- ncorrection	
    if(details) {
      attr(out$estimate, "details")$vars.used <- varname
	  attr(out$estimate, "details")$specimens.used <- specnames
    }      
    #class(out) <- c("dimorphEstDF", "data.frame")
    #class(out) <- "dimorphEstDF"
  }}
  return(out)
}


#' @noRd
dimorphMulti <- function(x, methodMulti="GMM", methodUni="SSD", sex=NULL, 
                         sex.female=1, center="geomean", ads=NULL,
                         templatevar=NULL, na.rm=T, ncorrection=F, details=F, dfout=F) {
  methodMulti <- match.arg(methodMulti, choices=c("GMM", "GMsize", "TM"))
  methodUni <- match.arg(methodUni, choices=c("SSD", "MMR", "BDI", "MoM", "FMA", "BFM",
                                              "ERM", "CV", "CVsex", "sdlog", "sdlogsex"))
  center <- match.arg(center, choices=c("geomean", "mean"))
  if (methodUni %in% c("CV", "CVsex")) center <- "mean"
  if (methodUni %in% c("sdlog", "sdlogsex")) center <- "geomean"
  propF.overall <- NA
  propF.realized <- NA
  if(!is.null(ads)) {
    x <- x[ads,]
	if(!is.null(sex)) sex <- sex[ads]
  }
  if (class(x)[1]=="data.frame") x <- as.matrix(x)
  if (class(x)[1]=="matrix") {
    if (inherits(x[,1], "character")) stop("All data in x must be numeric.")
    #if (class(x[,1])=="character") stop("All data in x must be numeric.")
  }
  {if (class(x)[1]=="list") {
    sex <- NULL
    propmissoverall <- NA
	nspecoverall <- NA
	nvarsoverall <- length(x)
  }
  else {
    if (is.null(rownames(x))) rownames(x) <- paste0("specimen_", 1:nrow(x))
    if (is.null(colnames(x))) colnames(x) <- paste0("VAR_", 1:ncol(x))
	nspecoverall <- nrow(x)
	nvarsoverall <- ncol(x)
    propmissoverall <- sum(is.na(x)) / (nrow(x)*ncol(x))
    # drop cases missing all measurements, and drop them from sex, too
    keep <- !(apply(is.na(x), prod, MARGIN=1))
    x <- x[keep,]
    if (!is.null(sex)) {
      sex <- sex[keep]
      if (!inherits(sex, "factor")) sex <- as.factor(sex)
      #if (!class(sex)=="factor") sex <- as.factor(sex)
      sex <- droplevels(sex)
      if (nlevels(sex)>2) stop("The sex variable must have exactly two levels.")
      if (!length(sex)==nrow(x)) stop("The number of rows in x and the length of the sex vector must be equal.")
      if (!(sex.female==1|sex.female==2)) stop("sex.female must be either 1 or 2.")
	  propF.overall <- sum(sex==levels(sex)[sex.female]) / length(sex)
    }
  }}
  { # start if loops for specific methodMulti
    if(methodMulti=="GMsize") { # start methodMulti "GMsize" (calculate geomean for each specimen and univariate dimorphism)
      if (class(x)[1]=="list") stop("'x' must be formatted as a matrix or data frame to use method 'GMsize'.")
      y <- apply(x, 1, dimorph::geomean, na.rm=F)
      if(sum(is.na(y)) > 0) warning("The geometric mean could not be calculated for all specimens.")
      if(sum(!is.na(y))==0) stop("No specimens are complete for all measurements.")
      if(sum(!is.na(y))< 2) stop("Fewer than two specimens are complete for all measurements.")
      out <- dimorph:::dimorphUni(y, method=methodUni, sex=sex, sex.female=sex.female, center=center,
                                 na.rm=na.rm, ncorrection=ncorrection, details=details, dfout=dfout)
	  #if(!is.null(sex)) propF.realized <- sum(sex[stats::complete.cases(y)]==levels(sex)[sex.female]) / sum(stats::complete.cases(y))
	  {if (!dfout) {
	    names(out) <- paste0(methodMulti, ".", methodUni)
        attr(out, "details")$n.vars.overall <- nvarsoverall
        attr(out, "details")$n.specimens.overall <- nspecoverall
	    attr(out, "details")$proportion.female.overall <- propF.overall
        attr(out, "details")$n.vars.realized <- ncol(x)
	    attr(out, "details")$proportion.missingdata.overall <- propmissoverall
	    attr(out, "details")$proportion.missingdata.realized <- 0
	    attr(out, "details")$proportion.templated <- NA
	    attr(out, "details")$methodMulti <- methodMulti
        if (details) {
          if(!is.null(colnames(x))) attr(out, "details")$vars.used <- colnames(x)
        } 
      }	  
      else {
	    rownames(out) <- paste0(methodMulti, ".", methodUni)
		out$methodMulti <- methodMulti
        out$n.vars.overall <- nvarsoverall
        out$n.specimens.overall <- nspecoverall
	    out$proportion.female.overall <- propF.overall
        out$n.vars.realized <- ncol(x)
		out$proportion.missingdata.overall <- propmissoverall
		out$proportion.missingdata.realized <- 0
		out$proportion.templated <- NA
        if (details) {
          if(!is.null(colnames(x))) attributes(out[[1]])$details$vars.used <- colnames(x)
        } 
      }}
      return(out)
    } # end methodMulti "GMsize"
    else if(methodMulti=="GMM") { # start methodMulti "GMM"
      #if (methodUni %in% c("CV", "CVsex", "sdlog", "sdlogsex")) warning("The 'GMM' multivariate method should NOT be used with the univariate methods 'CV', 'CVsex', 'sdlog', or 'sdlogsex'.")
      if (methodUni %in% c("CV", "CVsex", "sdlog", "sdlogsex")) stop("The 'GMM' multivariate method can NOT be used with the univariate methods 'CV', 'CVsex', 'sdlog', or 'sdlogsex'.")
      if (class(x)[1]=="matrix") x <- data.frame(x)
      #atleasttwo <- lapply(lapply(x, dimorph::not.is.na), sum) > 1
      atleasttwo <- lapply(lapply(x, function (y) return(!is.na(y))), sum) > 1
      if(sum(atleasttwo) < length(atleasttwo)) {
        warning(paste0("The following variable(s) were removed because they did not include\nat least two measurements:\n",
                       paste(names(atleasttwo)[!atleasttwo], collapse=", ")))
        x <- x[atleasttwo]
      }
	  # BFM needs at least three measurements per variable
	  if (methodUni=="BFM") {
        #atleast3 <- lapply(lapply(x, dimorph::not.is.na), sum) > 2
        atleast3 <- lapply(lapply(x, function (y) return(!is.na(y))), sum) > 2
        if(sum(atleast3) < length(atleast3)) {
          warning(paste0("The following variable(s) were removed because they did not include\nat least three measurements (required for BFM):\n", paste(names(atleast3)[!atleast3], collapse=", ")))
          x <- x[atleast3]
        }
	  }
	  #
	  # Drop specimens that no longer have any measurements
	  if (!"list" %in% class(x)) {
	    keep <- !(apply(is.na(x), prod, MARGIN=1))
	    if (sum(!keep) > 0) {
          x <- x[keep,]
          if(!is.null(sex)) sex <- sex[keep]
		  propF.realized <- sum(sex==levels(sex)[sex.female]) / length(sex)
          warning(paste0("The following specimen(s) were removed because they did not include\nat least one measurement after variables were removed:\n", paste(names(keep)[!keep], collapse=", ")))
	    }
	  }
	  if(!is.null(sex)) propF.realized <- sum(sex==levels(sex)[sex.female]) / length(sex)
      #
      ssd <- lapply(x, FUN=dimorph:::dimorphUni, method=methodUni, sex=sex,
                    sex.female=sex.female, center=center,
                    na.rm=na.rm, ncorrection=ncorrection, details=details, dfout=F)
      out <- dimorph::geomean(unlist(ssd), na.rm=F) # if some variables have NAs for estimates, no overall estimate should be returned
	  {if (!dfout) {
        attributes(out) <- attributes(ssd[[1]])
	    names(out) <- paste0(methodMulti, ".", methodUni)
	    attr(out, "details")$vars.used <- NA
	    attr(out, "details")$specimens.used <- NA
		{if ("list" %in% class(x)) {
          attr(out, "details")$n.vars.overall <- nvarsoverall
          attr(out, "details")$n.specimens.overall <- NA
          attr(out, "details")$n.vars.realized <- length(x)
          attr(out, "details")$n.specimens.realized <- NA
	      attr(out, "details")$proportion.missingdata.realized <- NA
          if (details) {
            attr(out, "details")$vars.used <- names(x)
          }
		}
		else {
          attr(out, "details")$n.vars.overall <- nvarsoverall
          attr(out, "details")$n.specimens.overall <- nspecoverall
          attr(out, "details")$n.vars.realized <- ncol(x)
          attr(out, "details")$n.specimens.realized <- nrow(x)
	      attr(out, "details")$proportion.missingdata.realized <- sum(is.na(x)) / (nrow(x)*ncol(x))
          if (details) {
            if(!is.null(colnames(x))) attr(out, "details")$vars.used <- colnames(x)
	        if(!is.null(rownames(x))) attr(out, "details")$specimens.used <- rownames(x)
          }
		}}
	    attr(out, "details")$proportion.missingdata.overall <- propmissoverall
	    attr(out, "details")$proportion.templated <- NA
	    attr(out, "details")$proportion.female.overall <- propF.overall
	    attr(out, "details")$proportion.female.realized <- propF.realized
	    attr(out, "details")$ratio.means <- NA
	    attr(out, "details")$methodMulti <- methodMulti
	    attr(out, "details")$model.parameters <- NA
		class(out) <- c("dimorphEst", "numeric")
      }	  
      else {
	    {if ("list" %in% class(x)) {
		  out <- data.frame(estimate=out,
			                methodUni=methodUni,
						    methodMulti=methodMulti,
							center=center,
							n.vars.overall=nvarsoverall,
							n.specimens.overall=nspecoverall,
							proportion.female.overall=propF.overall,
		                    n.vars.realized=length(x),
		                    n.specimens.realized=NA,
						    proportion.female.realized=propF.realized,
							proportion.missingdata.overall=propmissoverall,
						    proportion.missingdata.realized=NA,
						    proportion.templated=NA)
		}
		else {
		  out <- data.frame(estimate=out,
			                methodUni=methodUni,
						    methodMulti=methodMulti,
							center=center,
							n.vars.overall=nvarsoverall,
							n.specimens.overall=nspecoverall,
						    proportion.female.overall=propF.overall,
		                    n.vars.realized=ncol(x),
		                    n.specimens.realized=nrow(x),
						    proportion.female.realized=propF.realized,
						    proportion.missingdata.overall=propmissoverall,
						    proportion.missingdata.realized=sum(is.na(x)) / (nrow(x)*ncol(x)),
						    proportion.templated=NA)
		}}
		rownames(out) <- paste0(methodMulti, ".", methodUni)
        dets <- list(ratio.means=NA,
		             vars.used=NA,
					 specimens.used=NA,
					 model.parameters=NA,
					 estvalues="raw")	
		if (details) {
          if ("list" %in% class(x)) {
		    dets$vars.used <- names(x)
		  }
		  else {
		    if(!is.null(colnames(x))) dets$vars.used <- colnames(x)
            if(!is.null(rownames(x))) dets$specimens.used <- rownames(x)
		  }
        }
		attr(out[[1]], "details") <- dets
		class(out) <- c("dimorphEstDF", "data.frame")
		#class(out) <- "dimorphEstDF"
      }}
      return(out)
    } # end methodMulti "GMM"
	else if(methodMulti=="TM") { # start methodMulti "TM"
      if (class(x)[1]=="list") stop("'x' must be formatted as a matrix or data frame to use method 'TM'.")
      if (class(x)[1]=="matrix") x <- data.frame(x)
	  tmpEst <- dimorph:::TM(x, templatevar=templatevar)
      varsmissing <- colnames(tmpEst$template)[is.na(tmpEst$template)]
      {if(length(varsmissing)==ncol(tmpEst$templat)) {
	    warning("No specimen could act as template.")
	  }
      else if(length(varsmissing) > 0) {
        warning(paste0("The following variable(s) were removed because they\nwere not included in the template specimen:\n",
                       paste(varsmissing, collapse="\n")))
      }}
      specmissing <- names(tmpEst$TM)[is.na(tmpEst$TM)]
      if(length(specmissing) > 0) {
        warning(paste0("The following specimens(s) were removed because they\ndid not have templatable variables:\n",
                       paste(specmissing, collapse="\n")))
      }
      out <- dimorph:::dimorphUni(tmpEst$TM, method=methodUni, sex=sex,
                    sex.female=sex.female, center=center, ads=ads,
                    na.rm=na.rm, ncorrection=ncorrection, details=details, dfout=dfout)
	  x <- x[!is.na(tmpEst$template), !is.na(tmpEst$template)]
	  {if (!dfout) {
	    names(out) <- paste0(methodMulti, ".", methodUni)
        attr(out, "details")$n.vars.overall <- nvarsoverall
        attr(out, "details")$n.specimens.overall <- nspecoverall
		attr(out, "details")$proportion.female.overall <- propF.overall
        attr(out, "details")$n.vars.realized <- sum(!is.na(tmpEst$template))
	    attr(out, "details")$proportion.missingdata.overall <- propmissoverall
	    attr(out, "details")$proportion.missingdata.realized <- sum(is.na(x)) / (nrow(x)*ncol(x))
	    attr(out, "details")$proportion.templated <- tmpEst$prop.templated
	    attr(out, "details")$methodMulti <- methodMulti
	    attr(out, "details")$vars.used <- NA
	    attr(out, "details")$specimens.used <- NA
        if (details) {
          if(!is.null(colnames(x))) attr(out, "details")$vars.used <- colnames(x)
	      if(!is.null(rownames(x))) attr(out, "details")$specimens.used <- rownames(x)
        }
      }	  
      else {
	    rownames(out) <- paste0(methodMulti, ".", methodUni)
		out$methodMulti <- methodMulti
        out$n.vars.overall <- nvarsoverall
        out$n.specimens.overall <- nspecoverall
		out$proportion.female.overall <- propF.overall
        out$n.vars.realized <- sum(!is.na(tmpEst$template))
		out$proportion.missingdata.overall <- propmissoverall
		out$proportion.missingdata.realized <- sum(is.na(x)) / (nrow(x)*ncol(x))
		out$proportion.templated <- tmpEst$prop.templated
		attributes(out[[1]])$details$vars.used <- NA
		attributes(out[[1]])$details$specimens.used <- NA
        if (details) {
          if(!is.null(colnames(x))) attributes(out[[1]])$details$vars.used <- colnames(x)
          if(!is.null(rownames(x))) attributes(out[[1]])$details$specimens.used <- rownames(x)
        } 
      }}
      return(out)
    } # end methodMulti "TM"
  } # end if loops for specific methodMulti
}

#' @noRd
checktemplate <- function(template.ad, dat, miss) { # required for template method
  vars <- colnames(dat)[!is.na(dat[template.ad,])]
  n <- 0
  if (length(miss) > 1) n <- sum(apply(!is.na(dat[miss,vars]), 1, sum) > 0)
  else if (sum(!is.na(dat[miss,vars])) > 0) n <- 1
  return(n)
}

#' @noRd
TM <- function(x, templatevar=NULL) { # template method
  if (is.null(templatevar)) stop("'templatevar' must be specified with the name or\ncolumn number of a variable to estimate.")
  if (is.null(colnames(x))) colnames(x) <- paste0("VAR_", 1:ncol(x))
  if (inherits(templatevar, c("integer", "numeric"))) templatevar <- colnames(x)[templatevar]
  #if (class(templatevar)=="integer" | class(templatevar)=="numeric") templatevar <- colnames(x)[templatevar]
  x <- as.matrix(x)
  x.new <- x[,templatevar]
  n.templated <- 0
  NAout <- list(TM=x.new,
                prop.templated=0,
				template=x[1,,drop=F])
  NAout$template[1,] <- NA
  rownames(NAout$template) <- NULL
  if (sum(!is.na(x.new))==0) return(NAout)
  # get number of variables each specimen has, but zero if missing the template var
  counts <- apply(!(apply(x, 1, is.na)), 2, sum) * !is.na(x[,templatevar])
  # get adresses for specimens missing templatevar
  miss <- which(is.na(x.new))
  # Find specimen(s) that can template the most unknown specimens
  possibles <- which(counts > 1)
  if (length(possibles)==0) return(NAout)
  n.find <- apply(data.frame(possibles), 1, dimorph:::checktemplate, dat=x, miss=miss)
  template.ad <- possibles[which(n.find==max(n.find))]
  n.templated <- max(n.find)
  # If a tie for the most specimens, randomly choose one template specimen
  template.ad <- template.ad[sample(length(template.ad), 1)]
  template <- x[template.ad,]
  template <- template[!is.na(template)]
  for(zz in miss) {
    obs <- x[zz,]
	if(length(which(!is.na(obs)))==0) next
	mmt <- names(which(!is.na(obs)))
	mmt <- mmt[mmt %in% names(template)]
	if(length(mmt)==0) next
	mmt <- mmt[sample(length(mmt), 1)]
	x.new[zz] <- template[templatevar]/template[mmt] * x[zz, mmt]
  }
  prop.templated <- (sum(!is.na(x.new))-sum(!is.na(x[,templatevar])))/sum(!is.na(x.new))
  return(list(TM=x.new, prop.templated=prop.templated, template=x[template.ad,,drop=F]))
} # end function "TM"
