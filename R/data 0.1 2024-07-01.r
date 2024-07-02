#' Finite Mixture Analysis table from Godfrey et al. (1993)
#' @description Values of \emph{k} for a given sample size for the calculation of distribution
#'   means for a Gaussian finite mixture using equations in Godfrey et al. (1993). \emph{k}
#'   gives the mean expected number of standard deviations (as measured for the whole dataset)
#'   represented in the observed data range for a given sample size.
#' @docType data
#' @usage data(FMAtable)
#' @keywords datasets
#' @references Godfrey LR, Lyon SK, Sutherland MR. (1993) Sexual dimorphism in large-bodied 
#'   primates: The case of the subfossil lemurs. \emph{American Journal of Physical 
#'   Anthropology}. 90:315-334. (\href{https://doi.org/10.1002/ajpa.1330900306}{https://doi.org/10.1002/ajpa.1330900306})
#' @examples
#' data(FMAtable)
#' FMAtable
"FMAtable"


#' @title Long bone articular dimensions for four hominoid taxa
#' @description Postcranial metric data collected for western lowland gorillas, modern humans, common chimpanzees, and lar gibbons.  Further details can be found in Gordon (in review).
#' @format A data frame with 376 rows and 15 variables:
#' \describe{
#'   \item{\code{Species}}{\code{factor} Species name, with four possible levels (in decreasing order
#'         of postcranial dimorphism: \emph{\code{Gorilla gorilla}}, \emph{\code{Homo sapiens}}, 
#'         \emph{\code{Pan troglodytes}}, and \emph{\code{Hylobates lar}})}
#'   \item{\code{Museum}}{\code{factor} Collection housing each specimen}
#'   \item{\code{Collection.ID}}{\code{character} Unique specimen identifier within each collection}
#'   \item{\code{Sex}}{\code{factor} two possible levels: \code{F} and \code{M}}
#'   \item{\code{Wild}}{\code{factor} four possible levels: \code{Yes}, \code{Unknown}, \code{Died in captivity}, \code{No}}
#'   \item{\code{Mass.kg}}{\code{numeric} Body mass in kilograms when recorded in museum records.  Collections often record
#'      weight in pounds, in which case data were converted to mass in kilograms.}
#'   \item{\code{FHSI}}{\code{numeric} Femoral head superoinferior diameter: the maximum superoinferior diameter of 
#'         the femoral head.}
#'   \item{\code{TPML}}{\code{numeric} Tibial plateau mediolateral width: the maximum mediolateral width of the 
#'         articular surface of the tibial plateau.}
#'   \item{\code{TPMAP}}{\code{numeric} Tibial plateau medial condyle anteroposterior length: the maximum 
#'         anteroposterior length of the articular surface of the medial condyle of the tibial plateau.}
#'   \item{\code{TPLAP}}{\code{numeric} Tibial plateau lateral condyle anteroposterior length: the maximum 
#'         anteroposterior length of the articular surface of the lateral condyle of the tibial plateau.}
#'   \item{\code{HHMaj}}{\code{numeric} Humeral head major axis diameter: treating the articular surface of the 
#'         humeral head as a partial oblate spheroid, this is the length of the major axis passing through 
#'         the oblate spheroid.}
#'   \item{\code{HHMin}}{\code{numeric} Humeral head minor axis diameter: the maximum width of the articular surface 
#'         of the humeral head perpendicular to the major axis.}
#'   \item{\code{RHMaj}}{\code{numeric} Radial head major axis diameter: treating the radial head in proximal view 
#'         as an ellipse, this is the length of the major axis passing through the ellipse.}
#'   \item{\code{RHMin}}{\code{numeric} Radial head minor axis diameter: the maximum width of the radial head 
#'         perpendicular to the major axis.}
#'   \item{\code{RDAP}}{\code{numeric} Distal radius anteroposterior width: the distance between the anterior 
#'         and posterior extents of the boundary between the lunate and scaphoid facets of the distal 
#'         articular surface of the radius.}
#'   \item{\code{RDML}}{\code{numeric} Distal radius mediolateral breadth: the maximum width of the distal 
#'         articular surface of the radius when the medial point of this dimension is constrained to 
#'         the midpoint of the curve of the articulation with the distal ulna.}
#'}
#' @references Gordon AD (in review) Interpreting statistical significance in hominin dimorphism: Power and Type I error rates for resampling tests of univariate and missing-data multivariate size dimorphism estimation methods in the fossil record. \emph{Journal of Human Evolution}.
#' @examples
#' data(apelimbart)
#' plot(log10(HHMaj) ~ log10(FHSI), data=apelimbart, 
#'      pch=(21:24)[Species], bg=c(NA, "#00000040")[Sex])
#' legend("bottomright", legend=levels(apelimbart$Species), pch=21:24)
#' plot(log10(FHSI) ~ log10(Mass.kg), data=apelimbart, 
#'      pch=(21:24)[Species], bg=c(NA, "#00000040")[Sex])
#' legend("bottomright", legend=levels(apelimbart$Species), pch=21:24)
"apelimbart"

#' @title Long bone articular dimensions for two simulated fossil hominoid samples with missing data
#' @description Postcranial metric data collected for two "fauxil" species. \emph{Fauxil} sp. 1 individuals
#'   are extant gorillas from the \code{\link[dimorph]{apelimbart}} dataset with some data removed, and 
#'   \emph{Fauxil} sp. 2 individuals are modern humans from the \code{\link[dimorph]{apelimbart}} dataset 
#'   with some data removed.
#' @format A data frame with 35 rows and 13 variables:
#' \describe{
#'   \item{\code{Species}}{\code{factor} Species name, with four possible levels (in decreasing order
#'         of postcranial dimorphism: \code{\emph{Fauxil} sp. 1} and \code{\emph{Fauxil} sp. 2})}
#'   \item{\code{Museum}}{\code{factor} Collection housing each specimen}
#'   \item{\code{Collection.ID}}{\code{character} Unique specimen identifier within each collection}
#'   \item{\code{FHSI}}{\code{numeric} Femoral head superoinferior diameter: the maximum superoinferior diameter of 
#'         the femoral head.}
#'   \item{\code{TPML}}{\code{numeric} Tibial plateau mediolateral width: the maximum mediolateral width of the 
#'         articular surface of the tibial plateau.}
#'   \item{\code{TPMAP}}{\code{numeric} Tibial plateau medial condyle anteroposterior length: the maximum 
#'         anteroposterior length of the articular surface of the medial condyle of the tibial plateau.}
#'   \item{\code{TPLAP}}{\code{numeric} Tibial plateau lateral condyle anteroposterior length: the maximum 
#'         anteroposterior length of the articular surface of the lateral condyle of the tibial plateau.}
#'   \item{\code{HHMaj}}{\code{numeric} Humeral head major axis diameter: treating the articular surface of the 
#'         humeral head as a partial oblate spheroid, this is the length of the major axis passing through 
#'         the oblate spheroid.}
#'   \item{\code{HHMin}}{\code{numeric} Humeral head minor axis diameter: the maximum width of the articular surface 
#'         of the humeral head perpendicular to the major axis.}
#'   \item{\code{RHMaj}}{\code{numeric} Radial head major axis diameter: treating the radial head in proximal view 
#'         as an ellipse, this is the length of the major axis passing through the ellipse.}
#'   \item{\code{RHMin}}{\code{numeric} Radial head minor axis diameter: the maximum width of the radial head 
#'         perpendicular to the major axis.}
#'   \item{\code{RDAP}}{\code{numeric} Distal radius anteroposterior width: the distance between the anterior 
#'         and posterior extents of the boundary between the lunate and scaphoid facets of the distal 
#'         articular surface of the radius.}
#'   \item{\code{RDML}}{\code{numeric} Distal radius mediolateral breadth: the maximum width of the distal 
#'         articular surface of the radius when the medial point of this dimension is constrained to 
#'         the midpoint of the curve of the articulation with the distal ulna.}
#'}
#' @seealso \code{\link[dimorph]{apelimbart}}
#' @examples
#' data(fauxil)
#' fauxil
"fauxil"
