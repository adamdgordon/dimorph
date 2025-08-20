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
#' @description Postcranial metric data collected for western lowland gorillas, modern humans, common chimpanzees, and lar 
#'   gibbons from Gordon (2025a).
#' @format A data frame with 376 rows and 15 variables:
#' \describe{
#'   \item{\code{Species}}{\code{factor} Species name, with four possible levels (in decreasing order
#'         of postcranial dimorphism: \code{"Gorilla gorilla"}, \code{"Homo sapiens"}, 
#'         \code{"Pan troglodytes"}, and \code{"Hylobates lar"})}
#'   \item{\code{Museum}}{\code{factor} Collection housing each specimen}
#'   \item{\code{Collection.ID}}{\code{character} Unique specimen identifier within each collection}
#'   \item{\code{Sex}}{\code{factor} two possible levels: \code{"F"} and \code{"M"}}
#'   \item{\code{Wild}}{\code{factor} four possible levels: \code{"Yes"}, \code{"Unknown"}, \code{"Died in captivity"}, \code{"No"}}
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
#' @references Gordon AD. (2025a) Interpreting statistical significance in hominin dimorphism: Power and Type I error 
#'   rates for resampling tests of univariate and missing-data multivariate size dimorphism estimation methods in the 
#'   fossil record. \emph{Journal of Human Evolution}. 199:103630. 
#'   (\href{https://doi.org/10.1016/j.jhevol.2024.103630}{https://doi.org/10.1016/j.jhevol.2024.103630})
#' @examples
#' data(apelimbart)
#' ggplot2::ggplot(apelimbart, ggplot2::aes(x=log10(FHSI), y=log10(HHMaj), 
#'                                          color=Species, shape=Sex)) +
#'   ggplot2::geom_point(size=2) +
#'   ggplot2::scale_shape_manual(values = c(19,3))
#' ggplot2::ggplot(apelimbart, ggplot2::aes(x=log10(FHSI), y=log10(Mass.kg),
#'                                          color=Species, shape=Sex)) +
#'   ggplot2::geom_point(size=2) +
#'   ggplot2::scale_shape_manual(values = c(19,3))
"apelimbart"

#' @title Long bone articular dimensions for two simulated fossil hominoid samples with missing data
#' @description Postcranial metric data collected for two "fauxil" species. \emph{Fauxil} sp. 1 individuals
#'   are extant gorillas from the \code{\link[dimorph]{apelimbart}} dataset with some data removed, and 
#'   \emph{Fauxil} sp. 2 individuals are modern humans from the \code{\link[dimorph]{apelimbart}} dataset 
#'   with some data removed.  All data are from Gordon (2025a).
#' @format A data frame with 35 rows and 13 variables:
#' \describe{
#'   \item{\code{Species}}{\code{factor} Species name, with two possible levels (in decreasing order
#'         of postcranial dimorphism: \code{"Fauxil sp. 1"} and \code{"Fauxil sp. 2"})}
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
#' @references Gordon AD. (2025a) Interpreting statistical significance in hominin dimorphism: Power and Type I error 
#'   rates for resampling tests of univariate and missing-data multivariate size dimorphism estimation methods in the 
#'   fossil record. \emph{Journal of Human Evolution}. 199:103630. 
#'   (\href{https://doi.org/10.1016/j.jhevol.2024.103630}{https://doi.org/10.1016/j.jhevol.2024.103630})
#' @examples
#' data(fauxil)
#' fauxil
"fauxil"

#' @title Limb length and articular dimensions for six hominoid taxa
#' @description Linear postcranial metric data from Gordon et al. (2020) collected for 
#'   Sumatran orangutans, Bornean orangutans, western lowland gorillas, modern humans, common chimpanzees, 
#'   \emph{Australopithecus afarensis}, and \emph{A. africanus}.
#' @format A data frame with 216 rows and 15 variables:
#' \describe{
#'   \item{\code{Taxon}}{\code{factor} Taxon name, with six possible levels : \code{"Pongo"},
#'         \code{"Gorilla"}, \code{"Homo"}, \code{"Pan"}, \code{"A. afarensis"}, 
#'         and \code{"A. africanus"})}
#'   \item{\code{Species}}{\code{factor} Species name, with seven possible levels : \code{"Pongo abelii"},
#'         \code{"Pongo pygmaeus"}, \code{"Gorilla gorilla"}, \code{"Homo sapiens"}, 
#'         \code{"Pan troglodytes"}, \code{"A. afarensis"}, and \code{"A. africanus"})}
#'   \item{\code{Sex}}{\code{factor} Specimen sex, with three possible levels: \code{"F"}, \code{"M"}, and \code{"U"}}
#'   \item{\code{HUMHEAD}}{\code{numeric} Maximum anteroposterior (AP) diameter of the humeral head 
#'         taken perpendicular to the shaft axis.}
#'   \item{\code{ELBOW0.5}}{\code{numeric} Square root of the product of capitular height and articular width of the distal 
#'         humerus. Capitular height was taken from the anteroproximal border of capitulum to the distoposterior 
#'         border along the midline. Articular width was taken across the anterior aspect of the articular 
#'         surface from the lateral border of the capitulum to the medial edge of the articular surface.}
#'   \item{\code{RADTV}}{\code{numeric} Mediolateral (ML) diameter of the radial head.}
#'   \item{\code{FEMHEAD}}{\code{numeric} Maximum superoinferior (SI) diameter of the femoral head.}
#'   \item{\code{FEMSHAFT0.5}}{\code{numeric} Square root of the product of the AP and transverse diameters of the femoral shaft, 
#'         taken just inferior to the lesser trochanter.}
#'   \item{\code{DISTFEM0.5}}{\code{numeric} Square root of the product of the biepicondylar and shaft AP diameters of the distal femur.}
#'   \item{\code{PROXTIB0.5}}{\code{numeric} Square root of the product of the AP and transverse diameters of the proximal tibia. 
#'         The AP diameter was taken with one jaw of the calipers on the line connecting the posterior surfaces 
#'         of the medial and lateral condyles and the other jaw on the most distant point on the medial condyle. 
#'         Transverse diameter was the distance between the most lateral point on the lateral condyle to the most 
#'         medial point on the medial condyle (perpendicular to the AP diameter).}
#'   \item{\code{DISTTIB0.5}}{\code{numeric} Square root of the product of the AP and transverse diameters of the distal tibia. The 
#'         AP diameter is the distance between the most anterior and posterior points of the talar facet in the 
#'         AP plane. Transverse diameter is the distance between the midline of the medial malleolus and the 
#'         midline of the most medial point of the talar facet before the fibular facet begins.}
#'   \item{\code{Humerus.length}}{\code{integer} Humeral length measured to the nearest mm.}
#'   \item{\code{Radius.length}}{\code{integer} Radial length measured to the nearest mm.}
#'   \item{\code{Femur.length}}{\code{integer} Femoral length measured to the nearest mm.}
#'   \item{\code{Tibia.length}}{\code{integer} Tibial length measured to the nearest mm.}
#'}
#' @references Gordon AD, et al. (2020) Limb proportions and positional behavior: revisiting the theoretical 
#'   and empirical underpinnings for locomotor reconstruction in \emph{Australopithecus africanus}. In Zipfel B, 
#'   Richmond BG, and Ward CV, eds.: \emph{Hominid Postcranial Remains from Sterkfontein, South Africa, 1936-1995}. 
#'   Advances in Human Evolution Series. Oxford University Press. pp. 321-334.
#'   (\href{https://doi.org/10.1093/oso/9780197507667.003.0017}{Book Chapter}) (\href{https://doi.org/10.1093/oso/9780197507667.005.0003}{Appendix III}) (\href{https://doi.org/10.1093/oso/9780197507667.005.0004}{Appendix IV})
#' @examples
#' data(Gordonetal2020)
#' ggplot2::ggplot(Gordonetal2020, ggplot2::aes(x=log10(FEMHEAD), y=log10(HUMHEAD),
#'                                              color=Taxon, shape=Sex)) +
#'   ggplot2::geom_point(size=2) +
#'   ggplot2::scale_shape_manual(values = c(19,3, 17))
#' ggplot2::ggplot(Gordonetal2020, ggplot2::aes(x=log10(FEMHEAD), y=log10(RADTV),
#'                                              color=Taxon, shape=Sex)) +
#'   ggplot2::geom_point(size=2) +
#'   ggplot2::scale_shape_manual(values = c(19,3, 17))
"Gordonetal2020"

#' @title Articular dimensions for five hominoid taxa
#' @description Linear postcranial metric data from Gordon (2025b), some of which originally 
#'   appeared in Gordon \emph{et al.} (2020), collected for 
#'   western lowland gorillas (\emph{Gorilla gorilla gorilla}), modern humans, central chimpanzees 
#'   (\emph{Pan troglodytes troglodytes}), \emph{Australopithecus afarensis}, and \emph{A. africanus}.  
#'   Stratigraphic data and dates for \emph{A. afarensis} material from Campisano (2007).
#' @format A data frame with 216 rows and 14 variables:
#' \describe{
#'   \item{\code{Taxon}}{\code{factor} Taxon name, with five possible levels : 
#'         \code{"Gorilla"}, \code{"Homo"}, \code{"Pan"}, \code{"A. afarensis"}, 
#'         and \code{"A. africanus"}}
#'   \item{\code{Species}}{\code{factor} Species name, with seven possible levels : 
#'         \code{"Gorilla gorilla"}, \code{"Homo sapiens"}, 
#'         \code{"Pan troglodytes"}, \code{"A. afarensis"}, and \code{"A. africanus"}}
#'   \item{\code{Sex}}{\code{factor} Specimen sex, with three possible levels: \code{"F"}, \code{"M"}, and \code{"U"}}
#'   \item{\code{HUMHEAD}}{\code{numeric} Maximum anteroposterior (AP) diameter of the humeral head 
#'         taken perpendicular to the shaft axis.}
#'   \item{\code{ELBOW0.5}}{\code{numeric} Square root of the product of capitular height and articular width of the distal 
#'         humerus. Capitular height was taken from the anteroproximal border of capitulum to the distoposterior 
#'         border along the midline. Articular width was taken across the anterior aspect of the articular 
#'         surface from the lateral border of the capitulum to the medial edge of the articular surface.}
#'   \item{\code{RADTV}}{\code{numeric} Mediolateral (ML) diameter of the radial head.}
#'   \item{\code{FEMHEAD}}{\code{numeric} Maximum superoinferior (SI) diameter of the femoral head.}
#'   \item{\code{FEMSHAFT0.5}}{\code{numeric} Square root of the product of the AP and transverse diameters of the femoral shaft, 
#'         taken just inferior to the lesser trochanter.}
#'   \item{\code{DISTFEM0.5}}{\code{numeric} Square root of the product of the biepicondylar and shaft AP diameters of the distal femur.}
#'   \item{\code{PROXTIB0.5}}{\code{numeric} Square root of the product of the AP and transverse diameters of the proximal tibia. 
#'         The AP diameter was taken with one jaw of the calipers on the line connecting the posterior surfaces 
#'         of the medial and lateral condyles and the other jaw on the most distant point on the medial condyle. 
#'         Transverse diameter was the distance between the most lateral point on the lateral condyle to the most 
#'         medial point on the medial condyle (perpendicular to the AP diameter).}
#'   \item{\code{DISTTIB0.5}}{\code{numeric} Square root of the product of the AP and transverse diameters of the distal tibia. The 
#'         AP diameter is the distance between the most anterior and posterior points of the talar facet in the 
#'         AP plane. Transverse diameter is the distance between the midline of the medial malleolus and the 
#'         midline of the most medial point of the talar facet before the fibular facet begins.}
#'   \item{\code{Stratum}}{\code{factor} Submember of the Hadar formation that is the most-likely source of the fossil according
#'         to Campisano (2007).  Abbreviations: SH, Sidi Hakoma; DD, Denen Dora; KH, Hada Hadar.}
#'   \item{\code{Age.old}}{\code{numeric} Age (in millions of years) of the lower bound of the submemeber.}
#'   \item{\code{Age.young}}{\code{numeric} Age (in millions of years) of the upper bound of the submemeber.}
#'}
#' @references Campisano CJ. (2007) Tephrostratigraphy and hominin paleoenvironments of the Hadar Formation, 
#'   Afar Depression, Ethiopia (Ph.D.). Rutgers, The State University of New Jersey. 
#'   \href{https://www.proquest.com/docview/304805803}{https://www.proquest.com/docview/304805803}
#' @references Gordon AD. (2025b) Sexual size dimorphism in \emph{Australopithecus}: 
#'   postcranial dimorphism differs significantly among \emph{Australopithecus afarensis}, 
#'   \emph{A. africanus}, and modern humans despite low-power resampling analyses. \emph{American 
#'   Journal of Biological Anthropology}. 187:e70093. 
#'   \href{https://onlinelibrary.wiley.com/doi/10.1002/ajpa.70093}{https://onlinelibrary.wiley.com/doi/10.1002/ajpa.70093}
#' @references Gordon AD, et al. (2020) Limb proportions and positional behavior: revisiting the theoretical 
#'   and empirical underpinnings for locomotor reconstruction in \emph{Australopithecus africanus}. In Zipfel B, 
#'   Richmond BG, and Ward CV, eds.: \emph{Hominid Postcranial Remains from Sterkfontein, South Africa, 1936-1995}. 
#'   Advances in Human Evolution Series. Oxford University Press. pp. 321-334.
#'   (\href{https://doi.org/10.1093/oso/9780197507667.003.0017}{Book Chapter}) (\href{https://doi.org/10.1093/oso/9780197507667.005.0003}{Appendix III}) (\href{https://doi.org/10.1093/oso/9780197507667.005.0004}{Appendix IV})
#' @examples
#' data(GordonAJBA)
#' ggplot2::ggplot(GordonAJBA, ggplot2::aes(x=log10(FEMHEAD), y=log10(HUMHEAD),
#'                                          color=Taxon, shape=Sex)) +
#'   ggplot2::geom_point(size=2) +
#'   ggplot2::scale_shape_manual(values = c(19,3, 17))
#' ggplot2::ggplot(GordonAJBA, ggplot2::aes(x=log10(FEMHEAD), y=log10(RADTV),
#'                                          color=Taxon, shape=Sex)) +
#'   ggplot2::geom_point(size=2) +
#'   ggplot2::scale_shape_manual(values = c(19,3, 17))
"GordonAJBA"
