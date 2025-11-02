## ----include = FALSE---------------------------------------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup-------------------------------------------------------------------------------------------------------
library(dimorph)

## ----echo=F------------------------------------------------------------------------------------------------------
knitr::opts_knit$set(root.dir = "..")

## ----------------------------------------------------------------------------------------------------------------
str(apelimbart)

## ----fig.width=5, fig.height=4, echo=FALSE-----------------------------------------------------------------------
ggplot2::ggplot(apelimbart, ggplot2::aes(x=log10(FHSI), y=log10(HHMaj),
                                         color=Species, shape=Sex)) +
  ggplot2::geom_point(size=2, alpha = 0.5) +
  ggplot2::scale_shape_manual(values = c(19, 3, 17))

## ----------------------------------------------------------------------------------------------------------------
str(GordonAJBA)

## ----fig.width=5, fig.height=4, echo=FALSE-----------------------------------------------------------------------
ggplot2::ggplot(GordonAJBA, ggplot2::aes(x=log10(FEMHEAD), y=log10(HUMHEAD),
                                         color=Taxon, shape=Sex)) +
  ggplot2::geom_point(size=2, alpha = 0.5) +
  ggplot2::scale_shape_manual(values = c(19,3, 17))

## ----echo=FALSE--------------------------------------------------------------------------------------------------
displayvars <- c("Species", "Sex", "HUMHEAD", "ELBOW0.5", "RADTV", "FEMHEAD", 
                 "FEMSHAFT0.5", "DISTFEM0.5", "PROXTIB0.5", "DISTTIB0.5")
kableExtra::kable_styling(knitr::kable(
  head(GordonAJBA[GordonAJBA$Species %in% c("A. afarensis"),
             displayvars])),
  bootstrap_options = "condensed", font_size=10)

## ----------------------------------------------------------------------------------------------------------------
gorillas <- apelimbart[apelimbart$Species=="Gorilla gorilla",]

## ----------------------------------------------------------------------------------------------------------------
gorSSD  <- dimorph(x=gorillas$FHSI, method="SSD", sex=gorillas$Sex)
gorSSD

## ----------------------------------------------------------------------------------------------------------------
summary(gorSSD)

## ----------------------------------------------------------------------------------------------------------------
summary(dimorph(x=gorillas$FHSI, method="MMR"))
summary(dimorph(x=gorillas$FHSI, method="BDI"))
summary(dimorph(x=gorillas$FHSI, method="ERM"))

## ----------------------------------------------------------------------------------------------------------------
summary(dimorph(x=gorillas$FHSI, method="FMA"))
summary(dimorph(x=gorillas$FHSI, method="MoM"))
summary(dimorph(x=gorillas$FHSI, method="BFM"))

## ----------------------------------------------------------------------------------------------------------------
summary(dimorph(x=gorillas$FHSI, method="CV"))
summary(dimorph(x=gorillas$FHSI, method="CV", ncorrection=TRUE))

## ----------------------------------------------------------------------------------------------------------------
summary(dimorph(x=gorillas$FHSI, method="CVsex", sex=gorillas$Sex))
summary(dimorph(x=gorillas$FHSI, method="CVsex", sex=gorillas$Sex, ncorrection=TRUE))

## ----------------------------------------------------------------------------------------------------------------
summary(dimorph(x=gorillas$FHSI, method="sdlog"))
summary(dimorph(x=gorillas$FHSI, method="sdlogsex", sex=gorillas$Sex))

## ----------------------------------------------------------------------------------------------------------------
allmethods <- rbind(dimorph(x=gorillas$FHSI, method="SSD", dfout=TRUE, sex=gorillas$Sex),
                    dimorph(x=gorillas$FHSI, method="MMR", dfout=TRUE),
                    dimorph(x=gorillas$FHSI, method="BDI", dfout=TRUE),
                    dimorph(x=gorillas$FHSI, method="ERM", dfout=TRUE),
                    dimorph(x=gorillas$FHSI, method="FMA", dfout=TRUE),
                    dimorph(x=gorillas$FHSI, method="MoM", dfout=TRUE),
                    dimorph(x=gorillas$FHSI, method="BFM", dfout=TRUE),
                    dimorph(x=gorillas$FHSI, method="CV", dfout=TRUE),
                    dimorph(x=gorillas$FHSI, method="CVsex", dfout=TRUE, sex=gorillas$Sex),
                    dimorph(x=gorillas$FHSI, method="sdlog", dfout=TRUE),
                    dimorph(x=gorillas$FHSI, method="sdlogsex", dfout=TRUE, sex=gorillas$Sex))

## ----echo=F------------------------------------------------------------------------------------------------------
kableExtra::kable_styling(knitr::kable(allmethods),
  bootstrap_options = "condensed", font_size=10)

## ----------------------------------------------------------------------------------------------------------------
SSDvars_apelimbart <- c("FHSI", "TPML", "TPMAP", "TPLAP", "HHMaj",
                        "HHMin", "RHMaj", "RHMin", "RDAP", "RDML")
Gg_GMsize <- dimorph(x=gorillas[,SSDvars_apelimbart], method="SSD", methodMulti="GMsize",
                     sex=gorillas$Sex, details=TRUE)
summary(Gg_GMsize)

## ----------------------------------------------------------------------------------------------------------------
Gg_GMM <- dimorph(x=gorillas[,SSDvars_apelimbart], method="SSD", methodMulti="GMM",
                  sex=gorillas$Sex, details=TRUE)
summary(Gg_GMM)

## ----------------------------------------------------------------------------------------------------------------
gibbons <- apelimbart[apelimbart$Species=="Hylobates lar",]
bothmethods <- rbind(dimorph(x=gibbons[,SSDvars_apelimbart], method="SSD", methodMulti="GMsize",
                             sex=gibbons$Sex, dfout=TRUE),
                     dimorph(x=gibbons[,SSDvars_apelimbart], method="SSD", methodMulti="GMM",
                             sex=gibbons$Sex, dfout=TRUE))

## ----echo=F------------------------------------------------------------------------------------------------------
kableExtra::kable_styling(knitr::kable(bothmethods),
  bootstrap_options = "condensed", font_size=10)

## ----------------------------------------------------------------------------------------------------------------
SSDvars <- c("HUMHEAD", "ELBOW0.5", "RADTV", "FEMHEAD",    
             "FEMSHAFT0.5", "DISTFEM0.5", "PROXTIB0.5", "DISTTIB0.5")
Aafarensis <- GordonAJBA[GordonAJBA$Species=="A. afarensis", SSDvars]

## ----echo=F------------------------------------------------------------------------------------------------------
kableExtra::kable_styling(knitr::kable(Aafarensis),
  bootstrap_options = "condensed", font_size=10)

## ----------------------------------------------------------------------------------------------------------------
Aa_GMM <- dimorph(x=Aafarensis, method="MMR", methodMulti="GMM", details=TRUE)
Aa_GMM

## ----------------------------------------------------------------------------------------------------------------
summary(Aa_GMM, verbose=TRUE)

## ----------------------------------------------------------------------------------------------------------------
Aa_TM  <- dimorph(x=Aafarensis, method="MMR", methodMulti="TM", 
                  templatevar="FEMHEAD", details=TRUE)
Aa_TM

## ----------------------------------------------------------------------------------------------------------------
summary(Aa_TM, verbose=TRUE)

## ----------------------------------------------------------------------------------------------------------------
Aa_various <- rbind(dimorph(x=Aafarensis, method="MMR", methodMulti="GMM", dfout=TRUE),
                    suppressWarnings(dimorph(x=Aafarensis, method="MMR", methodMulti="TM",
                            dfout=TRUE, templatevar="FEMHEAD")),
                    dimorph(x=Aafarensis, method="BDI", methodMulti="GMM", dfout=TRUE),
                    suppressWarnings(dimorph(x=Aafarensis, method="BDI", methodMulti="TM",
                            dfout=TRUE, templatevar="FEMHEAD")),
                    dimorph(x=Aafarensis, method="ERM", methodMulti="GMM", dfout=TRUE),
                    suppressWarnings(dimorph(x=Aafarensis, method="ERM", methodMulti="TM",
                            dfout=TRUE, templatevar="FEMHEAD")),
                    dimorph(x=Aafarensis, method="FMA", methodMulti="GMM", dfout=TRUE),
                    suppressWarnings(dimorph(x=Aafarensis, method="FMA", methodMulti="TM",
                            dfout=TRUE, templatevar="FEMHEAD")),
                    dimorph(x=Aafarensis, method="MoM", methodMulti="GMM", dfout=TRUE),
                    suppressWarnings(dimorph(x=Aafarensis, method="MoM", methodMulti="TM",
                            dfout=TRUE, templatevar="FEMHEAD")),
                    dimorph(x=Aafarensis, method="BFM", methodMulti="GMM", dfout=TRUE),
                    suppressWarnings(dimorph(x=Aafarensis, method="BFM", methodMulti="TM",
                            dfout=TRUE, templatevar="FEMHEAD")),
                    suppressWarnings(dimorph(x=Aafarensis, method="CV", methodMulti="TM",
                            dfout=TRUE, templatevar="FEMHEAD")),
                    suppressWarnings(dimorph(x=Aafarensis, method="sdlog", methodMulti="TM",
                            dfout=TRUE, templatevar="FEMHEAD")))

## ----echo=F------------------------------------------------------------------------------------------------------
kableExtra::kable_styling(knitr::kable(Aa_various),
  bootstrap_options = "condensed", font_size=10)

## ----------------------------------------------------------------------------------------------------------------
gor <- apelimbart[apelimbart$Species=="Gorilla gorilla",]
hom <- apelimbart[apelimbart$Species=="Homo sapiens",]
pan <- apelimbart[apelimbart$Species=="Pan troglodytes",]
hyl <- apelimbart[apelimbart$Species=="Hylobates lar",]

## ----------------------------------------------------------------------------------------------------------------
nResample <- 1000

## ----------------------------------------------------------------------------------------------------------------
meths <- c("SSD", "MMR", "BDI", "ERM", "FMA", "MoM", "BFM", "CV", "CVsex", "sdlog", "sdlogsex")

## ----eval=FALSE--------------------------------------------------------------------------------------------------
#  # these lines of code are not actually being evaluated
#  bootdimorph(gor$FHSI)
#  bootdimorph(gor[,"FHSI"])
#  bootdimorph(gor[,"FHSI", drop=FALSE])

## ----------------------------------------------------------------------------------------------------------------
set.seed(5782) # not necessary, but generates reproducible example - see technical note above
bootsUgor <- bootdimorph(gor[,"FHSI", drop=FALSE], sex=gor$Sex, methsUni=meths, nResamp=nResample)
set.seed(5782) # not necessary, but generates reproducible example - see technical note above
bootsUhom <- bootdimorph(hom[,"FHSI", drop=FALSE], sex=hom$Sex, methsUni=meths, nResamp=nResample)
set.seed(5782) # not necessary, but generates reproducible example - see technical note above
bootsUpan <- bootdimorph(pan[,"FHSI", drop=FALSE], sex=pan$Sex, methsUni=meths, nResamp=nResample)
set.seed(5782) # not necessary, but generates reproducible example - see technical note above
bootsUhyl <- bootdimorph(hyl[,"FHSI", drop=FALSE], sex=hyl$Sex, methsUni=meths, nResamp=nResample)

## ----------------------------------------------------------------------------------------------------------------
bootsUgor

## ----------------------------------------------------------------------------------------------------------------
str(bootsUgor)

## ----------------------------------------------------------------------------------------------------------------
confint(bootsUgor)

## ----------------------------------------------------------------------------------------------------------------
confint(bootsUgor, type="bias")

## ----------------------------------------------------------------------------------------------------------------
confint(bootsUgor, conf.level=0.8, alternative="greater")

## ----fig.width=7, fig.height=7-----------------------------------------------------------------------------------
plot(bootsUgor)

## ----fig.width=7, fig.height=7-----------------------------------------------------------------------------------
plot(bootsUgor, exclude="FMA")

## ----fig.width=7, fig.height=7-----------------------------------------------------------------------------------
plot(bootsUhom)

## ----fig.width=7, fig.height=7-----------------------------------------------------------------------------------
plot(bootsUpan)

## ----fig.width=7, fig.height=7-----------------------------------------------------------------------------------
plot(bootsUhyl)

## ----fig.width=7, fig.height=7-----------------------------------------------------------------------------------
plot(bootsUgor, type="bias")

## ----fig.width=7, fig.height=7-----------------------------------------------------------------------------------
plot(bootsUhyl, type="bias")

## ----echo=FALSE--------------------------------------------------------------------------------------------------
set.seed(472)

## ----fig.width=7, fig.height=7-----------------------------------------------------------------------------------
bootsUgorUnbalanced <- bootdimorph(gor[c(1:47, 71:80), "FHSI", drop=FALSE],
                                   sex=gor$Sex[c(1:47, 71:80)], methsUni=meths, nResamp=nResample)
plot(bootsUgorUnbalanced, type="bias")

## ----------------------------------------------------------------------------------------------------------------
bootsMgor <- bootdimorph(gor[,SSDvars_apelimbart], sex=gor$Sex, methsUni=meths,
                         methsMulti=c("GMM", "GMsize"), nResamp=nResample)

## ----------------------------------------------------------------------------------------------------------------
bootsMgor

## ----fig.width=7, fig.height=7-----------------------------------------------------------------------------------
plot(bootsMgor)

## ----fig.width=7, fig.height=7-----------------------------------------------------------------------------------
plot(bootsMgor, exclude="FMA")

## ----fig.width=7, fig.height=7-----------------------------------------------------------------------------------
plot(bootsMgor, exclude="FMA", excludeMulti="GMsize")

## ----fig.width=7, fig.height=7-----------------------------------------------------------------------------------
plot(bootsMgor, type="bias")

## ----echo=F------------------------------------------------------------------------------------------------------
set.seed(32)

## ----fig.width=7, fig.height=7-----------------------------------------------------------------------------------
test_fossil_uni <- SSDtest(
     fossil=list("A. afarensis"=GordonAJBA[GordonAJBA$Species=="A. afarensis", "FEMHEAD", drop=FALSE]),
     comp=list("G. gorilla"=GordonAJBA[GordonAJBA$Species=="Gorilla gorilla", "FEMHEAD", drop=FALSE],
               "H. sapiens"=GordonAJBA[GordonAJBA$Species=="Homo sapiens", "FEMHEAD", drop=FALSE],
               "P. troglodytes"=GordonAJBA[GordonAJBA$Species=="Pan troglodytes", "FEMHEAD", drop=FALSE]),
     fossilsex=NULL,
     compsex=list("G. gorilla"=GordonAJBA[GordonAJBA$Species=="Gorilla gorilla", "Sex"],
                  "H. sapiens"=GordonAJBA[GordonAJBA$Species=="Homo sapiens", "Sex"],
                  "P. troglodytes"=GordonAJBA[GordonAJBA$Species=="Pan troglodytes", "Sex"]),
     methsUni=c("SSD", "MMR", "BDI"),
     nResamp=nResample)

## ----fig.width=7, fig.height=7-----------------------------------------------------------------------------------
test_fossil_uni

## ----------------------------------------------------------------------------------------------------------------
names(test_fossil_uni)

## ----------------------------------------------------------------------------------------------------------------
centers(test_fossil_uni)

## ----------------------------------------------------------------------------------------------------------------
pvals(test_fossil_uni)

## ----------------------------------------------------------------------------------------------------------------
suppressWarnings(pvals(test_fossil_uni, alternative="one.sided"))

## ----------------------------------------------------------------------------------------------------------------
names(test_fossil_uni$estimates)

## ----------------------------------------------------------------------------------------------------------------
test_fossil_uni$estimates[["A. afarensis"]]

## ----------------------------------------------------------------------------------------------------------------
test_fossil_uni$estimates[["G. gorilla"]]

## ----fig.width=7, fig.height=4-----------------------------------------------------------------------------------
plot(test_fossil_uni$estimates[["G. gorilla"]])

## ----fig.width=7, fig.height=4-----------------------------------------------------------------------------------
plot(test_fossil_uni)

## ----fig.width=7, fig.height=4-----------------------------------------------------------------------------------
plot(test_fossil_uni, est=2) # plots second method (MMR)

## ----fig.width=7, fig.height=4-----------------------------------------------------------------------------------
speciescolors <- c("A. afarensis"="#352A87", "A. africanus"="#F9FB0E", "G. gorilla"="#EABA4B",
                   "H. sapiens"="#09A9C0", "P. troglodytes"="#78BE7C")
plot(test_fossil_uni, est=2, groupcols=speciescolors) # change the colors

## ----fig.width=7, fig.height=3-----------------------------------------------------------------------------------
plot(test_fossil_uni, est=2, type="diff", # plots differences for all pairs of samples
     gridplot=FALSE, useradvance=FALSE) # prints each plot separately without waiting for user input

## ----fig.width=7, fig.height=3-----------------------------------------------------------------------------------
plot(test_fossil_uni, type="diff", est=2, diffs=c(1,2)) # plots diffs between first pair of samples

## ----echo=FALSE--------------------------------------------------------------------------------------------------
set.seed(8915)

## ----------------------------------------------------------------------------------------------------------------
test_2fossil_uni <- SSDtest(
     fossil=list("A. afarensis"=GordonAJBA[GordonAJBA$Species=="A. afarensis", "FEMHEAD"],
	               "A. africanus"=GordonAJBA[GordonAJBA$Species=="A. africanus", "FEMHEAD"]),
     comp=list("G. gorilla"=GordonAJBA[GordonAJBA$Species=="Gorilla gorilla", "FEMHEAD"],
               "H. sapiens"=GordonAJBA[GordonAJBA$Species=="Homo sapiens", "FEMHEAD"],
               "P. troglodytes"=GordonAJBA[GordonAJBA$Species=="Pan troglodytes", "FEMHEAD"]),
     fossilsex=NULL,
     compsex=list("G. gorilla"=GordonAJBA[GordonAJBA$Species=="Gorilla gorilla", "Sex"],
                  "H. sapiens"=GordonAJBA[GordonAJBA$Species=="Homo sapiens", "Sex"],
                  "P. troglodytes"=GordonAJBA[GordonAJBA$Species=="Pan troglodytes", "Sex"]),
     methsUni=c("MMR", "BDI"),
     exactcomp=FALSE,
     nResamp=nResample)

## ----fig.width=7, fig.height=7-----------------------------------------------------------------------------------
test_2fossil_uni

## ----fig.width=7, fig.height=4-----------------------------------------------------------------------------------
plot(test_2fossil_uni, est=1, groupcols=speciescolors)

## ----fig.width=7, fig.height=6-----------------------------------------------------------------------------------
plot(test_2fossil_uni, invert=2, est=1, groupcols=speciescolors)

## ----fig.width=7, fig.height=4-----------------------------------------------------------------------------------
plot(test_2fossil_uni, type="diff", est=2, diffs=c(2,3))

## ----fig.width=7, fig.height=7-----------------------------------------------------------------------------------
SSDvars <- c("HUMHEAD", "ELBOW0.5", "RADTV", "FEMHEAD",    
             "FEMSHAFT0.5", "DISTFEM0.5", "PROXTIB0.5", "DISTTIB0.5")

## ----echo=FALSE--------------------------------------------------------------------------------------------------
set.seed(5240)

## ----------------------------------------------------------------------------------------------------------------
test_fossil_multi <- SSDtest(
     fossil=list("A. afarensis"=GordonAJBA[GordonAJBA$Species=="A. afarensis", SSDvars]),
     comp=list("G. gorilla"=GordonAJBA[GordonAJBA$Species=="Gorilla gorilla", SSDvars],
               "H. sapiens"=GordonAJBA[GordonAJBA$Species=="Homo sapiens", SSDvars],
               "P. troglodytes"=GordonAJBA[GordonAJBA$Species=="Pan troglodytes", SSDvars]),
     fossilsex=NULL,
     compsex=list("G. gorilla"=GordonAJBA[GordonAJBA$Species=="Gorilla gorilla", "Sex"],
                  "H. sapiens"=GordonAJBA[GordonAJBA$Species=="Homo sapiens", "Sex"],
                  "P. troglodytes"=GordonAJBA[GordonAJBA$Species=="Pan troglodytes", "Sex"]),
     methsUni=c("MMR", "BDI"),
     methsMulti=c("GMM"),
     exactcomp = FALSE,
     datastruc="both",
     nResamp=nResample)

## ----------------------------------------------------------------------------------------------------------------
test_fossil_multi

## ----fig.width=7, fig.height=4-----------------------------------------------------------------------------------
plot(test_fossil_multi$estimates[["G. gorilla"]])

## ----fig.width=7, fig.height=4-----------------------------------------------------------------------------------
plot(test_fossil_multi)

## ----fig.width=7, fig.height=4-----------------------------------------------------------------------------------
plot(test_fossil_multi, est=2, # plot the 2nd method combination: MMR, GMM, "missing" datastructure
     groupcols=speciescolors) # use our previously specified color vector

## ----fig.width=7, fig.height=3-----------------------------------------------------------------------------------
plot(test_fossil_multi, est=2, type="diff",
     gridplot=FALSE, useradvance=FALSE)

## ----echo=FALSE--------------------------------------------------------------------------------------------------
set.seed(249)

## ----------------------------------------------------------------------------------------------------------------
test_2fossil_multi <- SSDtest(
     fossil=list("A. afarensis"=GordonAJBA[GordonAJBA$Species=="A. afarensis", SSDvars],
	             "A. africanus"=GordonAJBA[GordonAJBA$Species=="A. africanus", SSDvars]),
     comp=list("G. gorilla"=GordonAJBA[GordonAJBA$Species=="Gorilla gorilla", SSDvars],
               "H. sapiens"=GordonAJBA[GordonAJBA$Species=="Homo sapiens", SSDvars],
               "P. troglodytes"=GordonAJBA[GordonAJBA$Species=="Pan troglodytes", SSDvars]),
     fossilsex=NULL,
     compsex=list("G. gorilla"=GordonAJBA[GordonAJBA$Species=="Gorilla gorilla", "Sex"],
                  "H. sapiens"=GordonAJBA[GordonAJBA$Species=="Homo sapiens", "Sex"],
                  "P. troglodytes"=GordonAJBA[GordonAJBA$Species=="Pan troglodytes", "Sex"]),
     methsUni=c("MMR", "BDI"),
     methsMulti=c("GMM"),
     exactcomp = FALSE,
     exactfossil = TRUE,
     datastruc="both",
     nResamp=nResample)

## ----------------------------------------------------------------------------------------------------------------
test_2fossil_multi

## ----fig.width=7, fig.height=4-----------------------------------------------------------------------------------
plot(test_2fossil_multi, groupcols=speciescolors)

## ----fig.width=7, fig.height=6-----------------------------------------------------------------------------------
plot(test_2fossil_multi, est=2,
     groupcols=speciescolors, invert=c(1,2)) # invert first and second sample distributions

## ----fig.width=7, fig.height=3-----------------------------------------------------------------------------------
plot(test_2fossil_multi, est=2, type="diff",
     gridplot=FALSE, useradvance=FALSE)

## ----------------------------------------------------------------------------------------------------------------
citation(package="dimorph")

