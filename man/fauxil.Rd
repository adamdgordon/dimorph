% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_0.1_2025-07-16.r
\docType{data}
\name{fauxil}
\alias{fauxil}
\title{Long bone articular dimensions for two simulated fossil hominoid samples with missing data}
\format{
A data frame with 35 rows and 13 variables:
\describe{
  \item{\code{Species}}{\code{factor} Species name, with two possible levels (in decreasing order
        of postcranial dimorphism: \code{"Fauxil sp. 1"} and \code{"Fauxil sp. 2"})}
  \item{\code{Museum}}{\code{factor} Collection housing each specimen}
  \item{\code{Collection.ID}}{\code{character} Unique specimen identifier within each collection}
  \item{\code{FHSI}}{\code{numeric} Femoral head superoinferior diameter: the maximum superoinferior diameter of 
        the femoral head.}
  \item{\code{TPML}}{\code{numeric} Tibial plateau mediolateral width: the maximum mediolateral width of the 
        articular surface of the tibial plateau.}
  \item{\code{TPMAP}}{\code{numeric} Tibial plateau medial condyle anteroposterior length: the maximum 
        anteroposterior length of the articular surface of the medial condyle of the tibial plateau.}
  \item{\code{TPLAP}}{\code{numeric} Tibial plateau lateral condyle anteroposterior length: the maximum 
        anteroposterior length of the articular surface of the lateral condyle of the tibial plateau.}
  \item{\code{HHMaj}}{\code{numeric} Humeral head major axis diameter: treating the articular surface of the 
        humeral head as a partial oblate spheroid, this is the length of the major axis passing through 
        the oblate spheroid.}
  \item{\code{HHMin}}{\code{numeric} Humeral head minor axis diameter: the maximum width of the articular surface 
        of the humeral head perpendicular to the major axis.}
  \item{\code{RHMaj}}{\code{numeric} Radial head major axis diameter: treating the radial head in proximal view 
        as an ellipse, this is the length of the major axis passing through the ellipse.}
  \item{\code{RHMin}}{\code{numeric} Radial head minor axis diameter: the maximum width of the radial head 
        perpendicular to the major axis.}
  \item{\code{RDAP}}{\code{numeric} Distal radius anteroposterior width: the distance between the anterior 
        and posterior extents of the boundary between the lunate and scaphoid facets of the distal 
        articular surface of the radius.}
  \item{\code{RDML}}{\code{numeric} Distal radius mediolateral breadth: the maximum width of the distal 
        articular surface of the radius when the medial point of this dimension is constrained to 
        the midpoint of the curve of the articulation with the distal ulna.}
}
}
\usage{
fauxil
}
\description{
Postcranial metric data collected for two "fauxil" species. \emph{Fauxil} sp. 1 individuals
  are extant gorillas from the \code{\link[dimorph]{apelimbart}} dataset with some data removed, and 
  \emph{Fauxil} sp. 2 individuals are modern humans from the \code{\link[dimorph]{apelimbart}} dataset 
  with some data removed.  All data are from Gordon (2025a).
}
\examples{
data(fauxil)
fauxil
}
\references{
Gordon AD. (2025a) Interpreting statistical significance in hominin dimorphism: Power and Type I error 
  rates for resampling tests of univariate and missing-data multivariate size dimorphism estimation methods in the 
  fossil record. \emph{Journal of Human Evolution}. 199:103630. 
  (\href{https://doi.org/10.1016/j.jhevol.2024.103630}{https://doi.org/10.1016/j.jhevol.2024.103630})
}
\seealso{
\code{\link[dimorph]{apelimbart}}
}
\keyword{datasets}
