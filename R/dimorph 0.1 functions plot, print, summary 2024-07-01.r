#' @export
plot.dimorphResampledMulti <- function(x, plottype="draftsman", pch=21, pt.col=NULL,
                                       pt.bg="#00000030", na.cutoff=0.1, ...) {
  plottype <- match.arg(plottype, choices=c("draftsman", "sexratio"))
  tmp  <- x$estimates[,!colnames(x$estimates)%in%c("n", "nFem")]
  colrat="#0000FF"; colCV="#A52A2A"; colsd="#FFA500"
  if (plottype=="draftsman") {
    #ptmp <- apply(tmp, 2, sumNA)/nrow(tmp)
    ptmp <- apply(tmp, 2, function (x) return(sum(is.na(x))))/nrow(tmp)
    tmp  <- tmp[,(ptmp < na.cutoff)]
    plot(tmp, pch=pch, col=pt.col, bg=pt.bg, ...)
  }
  if (plottype=="sexratio") {
    if(is.na(x$estimate$nFem[1])) stop("Plot type 'sexratio' is only possible when sex data were included in the resampling procedure.")
    # Note that sex ratio can't be calculated for missing data when struc is a vector, so omit these
    if(x$sampleADS$strucClass=="vector") tmp <- tmp[,colnames(tmp)[-1*grep("missing", colnames(tmp))]]
    sxrat <- x$estimate$nFem/x$estimate$n
    varsrat <- colnames(tmp)
    varsCV  <- colnames(tmp)[grep("CV", colnames(tmp))]
    varssd  <- colnames(tmp)[grep("sdlog", colnames(tmp))]
    varsrat <- colnames(tmp)[!((varsrat %in% varsCV)|(varsrat %in% varssd))]
    nvars <- length(colnames(tmp))
    ngridcol <- ceiling(sqrt(nvars))
    ngridrow <- ceiling(nvars/ngridcol)
    parOLD <- par(no.readonly=T)
    par(mfrow=c(ngridrow, ngridcol), mar=c(2, 2, 2, 0) + 0.1)
    if (length(varsrat) > 0) {
      rngrat <- range(unlist(tmp[,varsrat]), na.rm=T)
      for (i in varsrat) {
        plot(x=sxrat, y=as.numeric(tmp[,i]), xlim=c(0,1), ylim=rngrat,
             pch=21, col=NULL, bg=paste0(colrat, "40"),
             main=paste0(i, "\nv. proportion of females"), xaxt="n",
             cex.main=0.9)
        axis(side=1, at=c(0,.25,0.5, 0.75, 1), labels=c("0.0", "", "0.5", "", "1.0"))
        #abline(h=0, col=lcol)
        #abline(h=1, col=lcol)
      }
    }
    if (length(varsCV) > 0) {
      rngCV <- range(unlist(tmp[,varsCV]), na.rm=T)
      for (i in varsCV) {
        plot(x=sxrat, y=as.numeric(tmp[,i]), xlim=c(0,1), ylim=rngCV,
             pch=21, col=NULL, bg=paste0(colCV, "40"),
             main=paste0(i, "\nv. proportion of females"), xaxt="n",
             cex.main=0.9)
        axis(side=1, at=c(0,.25,0.5, 0.75, 1), labels=c("0.0", "", "0.5", "", "1.0"))
        #abline(h=0, col=lcol)
        #abline(h=1, col=lcol)
      }
    }
    if (length(varssd) > 0) {
      rngsd <- range(unlist(tmp[,varssd]), na.rm=T)
      for (i in varssd) {
        plot(x=sxrat, y=as.numeric(tmp[,i]), xlim=c(0,1), ylim=rngsd,
             pch=21, col=NULL, bg=paste0(colsd, "40"),
             main=paste0(i, "\nv. proportion of females"), xaxt="n",
             cex.main=0.9)
        axis(side=1, at=c(0,.25,0.5, 0.75, 1), labels=c("0.0", "", "0.5", "", "1.0"))
        #abline(h=0, col=lcol)
        #abline(h=1, col=lcol)
      }
    }
    par(parOLD)
  } 
}


#' @export
hist.dimorphResampledMulti <- function(x, nbins=100, colrat="blue", colCV="brown", colsd="orange", ...) {
  tmp  <- x$estimates[,!colnames(x$estimates)%in%c("n", "nFem")]
  varsrat <- colnames(tmp)
  varsCV  <- colnames(tmp)[grep("CV", colnames(tmp))]
  varssd  <- colnames(tmp)[grep("sdlog", colnames(tmp))]
  varsrat <- colnames(tmp)[!((varsrat %in% varsCV)|(varsrat %in% varssd))]
  nvars <- length(colnames(tmp))
  ngridcol <- ceiling(sqrt(nvars))
  ngridrow <- ceiling(nvars/ngridcol)
  parOLD <- par(no.readonly=T)
  par(mfrow=c(ngridrow, ngridcol), mar=c(2, 1, 2, 0) + 0.1)
  lcol <- "#00000050"
  if (length(varsrat) > 0) {
    breaksrat <- hist(unlist(tmp[,varsrat]), breaks=nbins, plot=F)$breaks
    for (i in varsrat) {
      hist(as.numeric(tmp[,i]), breaks=breaksrat, col=colrat,
           main=paste0(i, "\nn valid resamples = ", length(sort(tmp[,i]))),
           yaxt="n", cex.main=0.9, border=F)
      abline(h=0, col=lcol)
    }
  }
  if (length(varsCV) > 0) {
    breaksCV <- hist(unlist(tmp[,varsCV]), breaks=nbins, plot=F)$breaks
    for (i in varsCV) {
      hist(as.numeric(tmp[,i]), breaks=breaksCV, col=colCV,
           main=paste0(i, "\nn valid resamples = ", length(sort(tmp[,i]))),
           yaxt="n", cex.main=0.9, border=F)
      abline(h=0, col=lcol)
    }
  }
  if (length(varssd) > 0) {
    breakssd <- hist(unlist(tmp[,varssd]), breaks=nbins, plot=F)$breaks
    for (i in varssd) {
      hist(as.numeric(tmp[,i]), breaks=breakssd, col=colsd,
           main=paste0(i, "\nn valid resamples = ", length(sort(tmp[,i]))),
           yaxt="n", cex.main=0.9, border=F)
      abline(h=0, col=lcol)
    }
  }
  par(parOLD) 
}

##' @export
print.dimorphEst <- function(x) {
  y <- x
  attributes(y) <- NULL
  names(y) <- names(x)
  if (attr(x,"details")$methodUni=="CV" | attr(x,"details")$methodUni=="CVsex") if(attr(x, "details")$ncorrection) names(y) <- paste0(names(y), "*")
  print(y)
}

##' @export
print.dimorphEstDF <- function(x) {
  attributes(x[[1]]) <- NULL
  class(x) <- "data.frame"
  print(x)
}

#' @export
summary.dimorphEst <- function(x, verbose=F) {
  txt <- paste0("estimate: ", round(x,5), "\n")
  txt <- paste0(txt, "univariate method: ", attr(x,"details")$methodUni)
  if (attr(x,"details")$methodUni=="CV" | attr(x,"details")$methodUni=="CVsex") if(attr(x, "details")$ncorrection)
        txt <- paste0(txt, " (sample size correction factor applied)")
  txt <- paste0(txt, "\n")
  if (!is.na(attr(x,"details")$methodMulti)) txt <- paste0(txt,"multivariate method: ",
	   attr(x,"details")$methodMulti,"\n")
  nvarsoveralltxt <- attr(x,"details")$n.vars.overall
  nvarsrealizedtxt <- attr(x,"details")$n.vars.realized
  if (!is.na(attr(x,"details")$methodMulti) & attr(x,"details")$methodMulti=="GMsize") {
    nvarsoveralltxt <- paste0("1 (geometric mean of ", nvarsoveralltxt, " variables)")
    nvarsrealizedtxt <- paste0("1 (geometric mean of ", nvarsrealizedtxt, " variables)")
  }
  propFoveralltxt <- round(attr(x,"details")$proportion.female.overall,5)
  propFrealizedtxt <- round(attr(x,"details")$proportion.female.realized,5)
  if(is.na(propFoveralltxt)) propFoveralltxt <- "unknown"
  if(is.na(propFrealizedtxt)) propFrealizedtxt <- "unknown"
  txt <- paste0(txt, "no. of variables (overall): ", nvarsoveralltxt,"\n")
  txt <- paste0(txt, "no. of specimens (overall): ", attr(x,"details")$n.specimens.overall,"\n")
  txt <- paste0(txt, "female proportion of sample (overall): ", propFoveralltxt,"\n")
  txt <- paste0(txt, "no. of variables (realized): ", nvarsrealizedtxt,"\n")
  txt <- paste0(txt, "no. of specimens (realized): ", attr(x,"details")$n.specimens.realized,"\n")
  txt <- paste0(txt, "female proportion of sample (realized): ", propFrealizedtxt,"\n")
  if(!is.na(attr(x,"details")$proportion.missingdata.overall))
    txt <- paste0(txt, "proportion of missing data (overall): ", round(attr(x,"details")$proportion.missingdata.overall,5),"\n")
  if(!is.na(attr(x,"details")$proportion.missingdata.realized))
    txt <- paste0(txt, "proportion of missing data (realized): ", round(attr(x,"details")$proportion.missingdata.realized,5),"\n")
  if(!is.na(attr(x,"details")$proportion.templated))
    txt <- paste0(txt, "proportion of template variable data estimated: ", round(attr(x,"details")$proportion.templated,5),"\n")
  centertxt <- attr(x,"details")$center
  if (!is.na(centertxt)) {
    if(centertxt=="geomean") centertxt <- "geometric mean"
    else if(centertxt=="mean") centertxt <- "arithmetic mean"
  }
  txt <- paste0(txt, "mean function: ", centertxt,"\n") 
  if (!is.na(attr(x,"details")$ratio.means[1])) txt <- paste0(txt,"ratio numerator and denominator: ",
	   paste(round(attr(x,"details")$ratio.means,2), collapse=", "),"\n")
  if (!is.na(attr(x,"details")$model.parameters[1])) {
    if (attr(x,"details")$methodUni=="BFM") {
      txt <- paste0(txt,"BFM model parameters:\n")
	  txt <- paste0(txt,"  BFM estimate of proportion of sample composed of smaller sex: ",
	         round(attr(x,"details")$model.parameters$pro[1],5),"\n")
	  modvar <- attr(x,"details")$model.parameters$variance$modelName
	  {if(modvar=="E") modvar <- "equal for both sexes"
	  else modvar <- "sex-specific"}
	  txt <- paste0(txt,"  BFM model of variance: ", modvar,"\n")
	  mvar <- round(attr(x,"details")$model.parameters$variance$sigmasq, 5)
	  if (attr(x,"details")$model.parameters$variance$modelName=="V") {
	    mvar <- paste0(mvar[1], " (smaller sex), ", mvar[2], " (larger sex)")
	  }
	  if (attr(x,"details")$center=="geomean") mvar <- paste0(mvar, " (logged data)")
	  if (attr(x,"details")$center=="mean") mvar <- paste0(mvar, " (raw data)")
	  txt <- paste0(txt,"  BFM estimate of variance: ", mvar,"\n")
	}
  }
  if(verbose) {
    if (!is.na(attr(x,"details")$vars.used[1])) txt <- paste0(txt, "\nIncluded variables:\n", 
                  paste(attr(x,"details")$vars.used, collapse="\n"), "\n")
    if (!is.na(attr(x,"details")$specimens.used[1])) txt <- paste0(txt, "\nIncluded specimens:\n", 
                  paste(attr(x,"details")$specimens.used, collapse="\n"), "\n")
  }
  cat(txt)
}

#' @export
summary.dimorphEstDF <- function(x, verbose=F) {
  txt <- paste0("estimate: ", round(x$estimate,5), "\n")
  txt <- paste0(txt, "univariate method: ", x$methodUni,"\n")
  if (!is.na(x$methodMulti)) txt <- paste0(txt,"multivariate method: ",
	   x$methodMulti,"\n")
  nvarsoveralltxt <- x$n.vars.overall
  nvarsrealizedtxt <- x$n.vars.realized
  if (!is.na(x$methodMulti) & x$methodMulti=="GMsize") {
    nvarsoveralltxt <- paste0("1 (geometric mean of ", nvarsoveralltxt, " variables)")
    nvarsrealizedtxt <- paste0("1 (geometric mean of ", nvarsrealizedtxt, " variables)")
  }
  propFoveralltxt <- round(x$proportion.female.overall,5)
  propFrealizedtxt <- round(x$proportion.female.realized,5)
  if(is.na(propFoveralltxt)) propFoveralltxt <- "unknown"
  if(is.na(propFrealizedtxt)) propFrealizedtxt <- "unknown"
  txt <- paste0(txt, "no. of variables (overall): ", nvarsoveralltxt,"\n")
  txt <- paste0(txt, "no. of specimens (overall): ", x$n.specimens.overall,"\n")
  txt <- paste0(txt, "female proportion of sample (overall): ", propFoveralltxt,"\n")
  txt <- paste0(txt, "no. of variables (realized): ", nvarsrealizedtxt,"\n")
  txt <- paste0(txt, "no. of specimens (realized): ", x$n.specimens.realized,"\n")
  txt <- paste0(txt, "female proportion of sample (realized): ", propFrealizedtxt,"\n")
  if(!is.na(x$proportion.missingdata.overall))
    txt <- paste0(txt, "proportion of missing data (overall): ", round(x$proportion.missingdata.overall,5),"\n")
  if(!is.na(x$proportion.missingdata.realized))
    txt <- paste0(txt, "proportion of missing data (realized): ", round(x$proportion.missingdata.realized,5),"\n")
  if(!is.na(x$proportion.templated))
    txt <- paste0(txt, "proportion of template variable data estimated: ", round(x$proportion.templated,5),"\n")
  centertxt <- x$center
  if (!is.na(centertxt)) {
    if(centertxt=="geomean") centertxt <- "geometric mean"
    else if(centertxt=="mean") centertxt <- "arithmetic mean"
  }
  txt <- paste0(txt, "mean function: ", centertxt,"\n") 
  if (!is.na(attr(x[[1]],"details")$ratio.means[1])) txt <- paste0(txt,"ratio numerator and denominator: ",
	   paste(round(attr(x[[1]],"details")$ratio.means,2), collapse=", "),"\n")
  if (!is.na(attr(x[[1]],"details")$model.parameters[1])) {
    if (x$methodUni=="BFM") {
      txt <- paste0(txt,"BFM model parameters:\n")
	  txt <- paste0(txt,"  BFM estimate of proportion of sample composed of smaller sex: ",
	         round(attr(x[[1]],"details")$model.parameters$pro[1],5),"\n")
	  modvar <- attr(x[[1]],"details")$model.parameters$variance$modelName
	  {if(modvar=="E") modvar <- "equal for both sexes"
	  else modvar <- "sex-specific"}
	  txt <- paste0(txt,"  BFM model of variance: ", modvar,"\n")
	  mvar <- round(attr(x[[1]],"details")$model.parameters$variance$sigmasq, 5)
	  if (attr(x[[1]],"details")$model.parameters$variance$modelName=="V") {
	    mvar <- paste0(mvar[1], " (smaller sex), ", mvar[2], " (larger sex)")
	  }
	  if (attr(x[[1]],"details")$center=="geomean") mvar <- paste0(mvar, " (logged data)")
	  if (attr(x[[1]],"details")$center=="mean") mvar <- paste0(mvar, " (raw data)")
	  txt <- paste0(txt,"  BFM estimate of variance: ", mvar,"\n")
	}
  }
  if(verbose) {
    if (!is.na(attr(x[[1]],"details")$vars.used[1])) txt <- paste0(txt, "\nIncluded variables:\n", 
                  paste(attr(x[[1]],"details")$vars.used, collapse="\n"), "\n")
    if (!is.na(attr(x[[1]],"details")$specimens.used[1])) txt <- paste0(txt, "\nIncluded specimens:\n", 
                  paste(attr(x[[1]],"details")$specimens.used, collapse="\n"), "\n")
  }
  cat(txt)
}


#' @export
plot.dimorphResampledUni <- function(x, type="estimate") {
  CIheight <- 0.4
  type <- match.arg(type, choices=c("estimate", "bias"))
  if (type=="bias" & sum(is.na(x$estimates$bias))==length(x$estimates$bias)) stop("No bias data are present.")
  x$estimates$method <- droplevels(x$estimates$methodUni:x$estimates$center)
  if (!is.null(x$CI)) x$CI$method <- droplevels(x$CI$methodUni:x$CI$center)
  if (!is.null(x$CIbias)) x$CIbias$method <- droplevels(x$CIbias$methodUni:x$CIbias$center)
  {if (type=="estimate") {
    SSDvars <- c("SSD", "MMR", "BDI", "ERM", "FMA", "MoM", "BFM")
    CVvars <- c("CV", "CVsex")
    sdlogvars <- c("sdlog", "sdlogsex")
  }
  else if (type=="bias") {
    SSDvars <- c("MMR", "BDI", "ERM", "FMA", "MoM", "BFM")
    CVvars <- c("CV")
    sdlogvars <- c("sdlog")
  }}
  nSSDv <- sum(unique(x$estimates$methodUni) %in% SSDvars)
  nCVv <- sum(unique(x$estimates$methodUni) %in% CVvars)
  nsdlogv <- sum(unique(x$estimates$methodUni) %in% sdlogvars)
  pltSSD <- NULL
  pltCV <- NULL
  pltsdlog <- NULL
  # generate pltSSD
  if (nSSDv > 0 & sum(is.na(x$estimates[x$estimates$methodUni %in% SSDvars,type])) < length(x$estimates[x$estimates$methodUni %in% SSDvars,type])) {
    {if (type=="estimate") {
	  pltSSD <- ggplot2::ggplot()
	  {if (attr(x$estimates, "estvalues")=="raw") pltSSD <- pltSSD + ggplot2::scale_x_log10() + ggplot2::geom_vline(xintercept=1) + ggplot2::xlab("estimate (axis scale is log10 transformed)")
	  else pltSSD <- pltSSD + ggplot2::geom_vline(xintercept=0)}
	  pltSSD <- pltSSD + ggplot2::geom_violin(data=x$estimates[x$estimates$methodUni %in% SSDvars,],
                          ggplot2::aes(y=method, x=estimate, fill=methodUni)) +
            ggplot2::stat_summary(data=x$estimates[x$estimates$methodUni %in% SSDvars,],
                          ggplot2::aes(y=method, x=estimate, fill=methodUni),fun=median, geom="point", size=2, color="black") +
            ggplot2::theme(legend.position = "none")
    }
	else if (type=="bias") pltSSD <- ggplot2::ggplot() + ggplot2::geom_vline(xintercept=0) +
            ggplot2::geom_violin(data=x$estimates[x$estimates$methodUni %in% SSDvars,],
                          ggplot2::aes(y=method, x=bias, fill=methodUni)) +
            ggplot2::stat_summary(data=x$estimates[x$estimates$methodUni %in% SSDvars,],
                          ggplot2::aes(y=method, x=bias, fill=methodUni),fun=median, geom="point", size=2, color="black") +
            ggplot2::theme(legend.position = "none") +
			ggplot2::xlab("bias from sample SSD")}
    CIdat <- NULL
    if (type=="estimate" & !is.null(x$CI)) CIdat <- x$CI		  
    else if (type=="bias" & !is.null(x$CIbias)) CIdat <- x$CIbias			  
    if (!is.null(CIdat)) {		  
      # only do one or the other if upper or lower bound
	  segdatSSDlow <- NULL
	  segdatSSDhigh <- NULL
	  if (attr(x, "alternative")=="two.sided" | attr(x, "alternative")=="greater") {
        segdatSSDlow <- data.frame(x=CIdat$lower_lim[CIdat$methodUni %in% SSDvars],
                                 y=CIheight,
								 xend=CIdat$lower_lim[CIdat$methodUni %in% SSDvars],
								 yend=-CIheight)
        segdatSSDlow$y <- segdatSSDlow$y + as.integer(rownames(segdatSSDlow))
        segdatSSDlow$yend <- segdatSSDlow$yend + as.integer(rownames(segdatSSDlow))
	  }
	  if (attr(x, "alternative")=="two.sided" | attr(x, "alternative")=="less") {
        segdatSSDhigh <- data.frame(x=CIdat$upper_lim[CIdat$methodUni %in% SSDvars],
                                 y=CIheight,
								 xend=CIdat$upper_lim[CIdat$methodUni %in% SSDvars],
								 yend=-CIheight)
        segdatSSDhigh$y <- segdatSSDhigh$y + as.integer(rownames(segdatSSDhigh))
        segdatSSDhigh$yend <- segdatSSDhigh$yend + as.integer(rownames(segdatSSDhigh))
	  }
	  segdatSSD <- rbind(segdatSSDlow, segdatSSDhigh)
      pltSSD <- pltSSD + ggplot2::geom_segment(data=segdatSSD, ggplot2::aes(x=x, y=y, xend=xend, yend=yend),
                                       color="black")
      pltSSD <- pltSSD + ggplot2::geom_segment(data=segdatSSD, ggplot2::aes(x=x, y=y, xend=xend, yend=yend),
                                       color="red", linetype=3)
    }
  } # end pltSSD
  # generate pltCV
  if (nCVv > 0 & sum(is.na(x$estimates[x$estimates$methodUni %in% CVvars,type])) < length(x$estimates[x$estimates$methodUni %in% CVvars,type])) {
    {if (type=="estimate") pltCV <- ggplot2::ggplot() + ggplot2::geom_vline(xintercept=0) +
            ggplot2::geom_violin(data=x$estimates[x$estimates$methodUni %in% CVvars,],
                          ggplot2::aes(y=method, x=estimate, fill=methodUni)) +
            ggplot2::stat_summary(data=x$estimates[x$estimates$methodUni %in% CVvars,],
                          ggplot2::aes(y=method, x=estimate, fill=methodUni),fun=median, geom="point", size=2, color="black") +
            ggplot2::theme(legend.position = "none")
	else if (type=="bias") pltCV <- ggplot2::ggplot() + ggplot2::geom_vline(xintercept=0) +
            ggplot2::geom_violin(data=x$estimates[x$estimates$methodUni %in% CVvars,],
                          ggplot2::aes(y=method, x=bias, fill=methodUni)) +
            ggplot2::stat_summary(data=x$estimates[x$estimates$methodUni %in% CVvars,],
                          ggplot2::aes(y=method, x=bias, fill=methodUni),fun=median, geom="point", size=2, color="black") +
            ggplot2::theme(legend.position = "none") +
			ggplot2::xlab("bias from sample CVsex")}
    CIdat <- NULL
    if (type=="estimate" & !is.null(x$CI)) CIdat <- x$CI		  
    else if (type=="bias" & !is.null(x$CIbias)) CIdat <- x$CIbias			  
    if (!is.null(CIdat)) {		  
      # only do one or the other if upper or lower bound
	  segdatCVlow <- NULL
	  segdatCVhigh <- NULL
	  if (attr(x, "alternative")=="two.sided" | attr(x, "alternative")=="greater") {
        segdatCVlow <- data.frame(x=CIdat$lower_lim[CIdat$methodUni %in% CVvars],
                                 y=CIheight,
								 xend=CIdat$lower_lim[CIdat$methodUni %in% CVvars],
								 yend=-CIheight)
        segdatCVlow$y <- segdatCVlow$y + as.integer(rownames(segdatCVlow))
        segdatCVlow$yend <- segdatCVlow$yend + as.integer(rownames(segdatCVlow))
	  }
	  if (attr(x, "alternative")=="two.sided" | attr(x, "alternative")=="less") {
        segdatCVhigh <- data.frame(x=CIdat$upper_lim[CIdat$methodUni %in% CVvars],
                                 y=CIheight,
								 xend=CIdat$upper_lim[CIdat$methodUni %in% CVvars],
								 yend=-CIheight)
        segdatCVhigh$y <- segdatCVhigh$y + as.integer(rownames(segdatCVhigh))
        segdatCVhigh$yend <- segdatCVhigh$yend + as.integer(rownames(segdatCVhigh))
	  }
	  segdatCV <- rbind(segdatCVlow, segdatCVhigh)
      pltCV <- pltCV + ggplot2::geom_segment(data=segdatCV, ggplot2::aes(x=x, y=y, xend=xend, yend=yend),
                                       color="black")
      pltCV <- pltCV + ggplot2::geom_segment(data=segdatCV, ggplot2::aes(x=x, y=y, xend=xend, yend=yend),
                                       color="red", linetype=3)
    }
  } # end pltCV
  # generate pltsdlog
  if (nsdlogv > 0 & sum(is.na(x$estimates[x$estimates$methodUni %in% sdlogvars,type])) < length(x$estimates[x$estimates$methodUni %in% sdlogvars,type])) {
    {if (type=="estimate") pltsdlog <- ggplot2::ggplot() + ggplot2::geom_vline(xintercept=0) +
            ggplot2::geom_violin(data=x$estimates[x$estimates$methodUni %in% sdlogvars,],
                          ggplot2::aes(y=method, x=estimate, fill=methodUni)) +
            ggplot2::stat_summary(data=x$estimates[x$estimates$methodUni %in% sdlogvars,],
                          ggplot2::aes(y=method, x=estimate, fill=methodUni),fun=median, geom="point", size=2, color="black") +
            ggplot2::theme(legend.position = "none")
    else if (type=="bias") pltsdlog <- ggplot2::ggplot() + ggplot2::geom_vline(xintercept=0) +
            ggplot2::geom_violin(data=x$estimates[x$estimates$methodUni %in% sdlogvars,],
                          ggplot2::aes(y=method, x=bias, fill=methodUni)) +
            ggplot2::stat_summary(data=x$estimates[x$estimates$methodUni %in% sdlogvars,],
                          ggplot2::aes(y=method, x=bias, fill=methodUni),fun=median, geom="point", size=2, color="black") +
            ggplot2::theme(legend.position = "none") +
			ggplot2::xlab("bias from sample sdlogsex")}
    CIdat <- NULL
    if (type=="estimate" & !is.null(x$CI)) CIdat <- x$CI		  
    else if (type=="bias" & !is.null(x$CIbias)) CIdat <- x$CIbias			  
    if (!is.null(CIdat)) {		  
      # only do one or the other if upper or lower bound
	  segdatsdloglow <- NULL
	  segdatsdloghigh <- NULL
	  if (attr(x, "alternative")=="two.sided" | attr(x, "alternative")=="greater") {
        segdatsdloglow <- data.frame(x=CIdat$lower_lim[CIdat$methodUni %in% sdlogvars],
                                 y=CIheight,
								 xend=CIdat$lower_lim[CIdat$methodUni %in% sdlogvars],
								 yend=-CIheight)
        segdatsdloglow$y <- segdatsdloglow$y + as.integer(rownames(segdatsdloglow))
        segdatsdloglow$yend <- segdatsdloglow$yend + as.integer(rownames(segdatsdloglow))
	  }
	  if (attr(x, "alternative")=="two.sided" | attr(x, "alternative")=="less") {
        segdatsdloghigh <- data.frame(x=CIdat$upper_lim[CIdat$methodUni %in% sdlogvars],
                                 y=CIheight,
								 xend=CIdat$upper_lim[CIdat$methodUni %in% sdlogvars],
								 yend=-CIheight)
        segdatsdloghigh$y <- segdatsdloghigh$y + as.integer(rownames(segdatsdloghigh))
        segdatsdloghigh$yend <- segdatsdloghigh$yend + as.integer(rownames(segdatsdloghigh))
	  }
	  segdatsdlog <- rbind(segdatsdloglow, segdatsdloghigh)
      pltsdlog <- pltsdlog + ggplot2::geom_segment(data=segdatsdlog, ggplot2::aes(x=x, y=y, xend=xend, yend=yend),
                                       color="black")
      pltsdlog <- pltsdlog + ggplot2::geom_segment(data=segdatsdlog, ggplot2::aes(x=x, y=y, xend=xend, yend=yend),
                                       color="red", linetype=3)
    }
  } # end pltsdlog
  # build plotting text
  n <- nrow(x$sampleADS$comparative)
  nResamp <- x$sampleADS$nResamp
  toptxt <- "Resampled"
  if (!is.null(attr(x, "resampling"))) {
    if (attr(x, "resampling")=="bootstrap") toptxt <- "Bootstrapped"
  }
  {if (type=="estimate") toptxt <- paste0(toptxt, " estimates")
  else if (type=="bias") toptxt <- paste0(toptxt, " estimate bias")}
  if (!is.null(attr(x, "conf.level"))) {
    toptxt <- paste0(toptxt, " with ")
	if (attr(x, "alternative")=="two.sided") toptxt <- paste0(toptxt, "two-sided ",
	     attr(x, "conf.level")*100, "% confidence intervals\n(univariate data, n=", n, ", n resampled datasets=", nResamp, ")")
	else if (attr(x, "alternative")=="less") toptxt <- paste0(toptxt, "one-sided ",
	     attr(x, "conf.level")*100, "% confidence intervals\n(upper bound, univariate data, n=", n, ", n resampled datasets=", nResamp, ")")
	else if (attr(x, "alternative")=="greater") toptxt <- paste0(toptxt, "one-sided ", 
	     attr(x, "conf.level")*100, "% confidence intervals\n(lower bound, univariate data, n=", n, ", n resampled datasets=", nResamp, ")")
  }
  {if (is.null(pltCV) & is.null(pltsdlog)) {
    return(gridExtra::grid.arrange(pltSSD, top=toptxt))
  }
  else if (is.null(pltSSD) & is.null(pltsdlog)) {
    return(gridExtra::grid.arrange(pltCV, top=toptxt))
  }
  else if (is.null(pltSSD) & is.null(pltCV)) {
    return(gridExtra::grid.arrange(pltsdlog, top=toptxt))
  }
  else if (is.null(pltSSD)) {
    return(gridExtra::grid.arrange(pltCV, pltsdlog,
	         heights=c(nCVv,nsdlogv), ncol = 1,
			 top=toptxt))
  }
  else if (is.null(pltCV)) {
    return(gridExtra::grid.arrange(pltSSD, pltsdlog,
	         heights=c(nSSDv,nsdlogv), ncol = 1,
			 top=toptxt))
  }
  else if (is.null(pltsdlog)) {
    return(gridExtra::grid.arrange(pltSSD, pltCV,
	         heights=c(nSSDv,nCVv), ncol = 1,
			 top=toptxt))
  }
  else {
    return(gridExtra::grid.arrange(pltSSD, pltCV, pltsdlog,
	         heights=c(nSSDv,nCVv,nsdlogv), ncol = 1,
			 top=toptxt))
  }}
}

#' @export
plot.dimorphResampledMulti <- function(x, type="estimate") {
  CIheight <- 0.4
  type <- match.arg(type, choices=c("estimate", "bias"))
  if (type=="bias" & sum(is.na(x$estimates$bias))==length(x$estimates$bias)) stop("No bias data are present.")
  x$estimates$method <- droplevels(x$estimates$methodUni:x$estimates$methodMulti:x$estimates$center:x$estimates$datastructure)
  if (!is.null(x$CI)) x$CI$method <- droplevels(x$CI$methodUni:x$CI$methodMulti:x$CI$center:x$CI$datastructure)
  if (!is.null(x$CIbias)) x$CIbias$method <- droplevels(x$CIbias$methodUni:x$CIbias$methodMulti:x$CIbias$center:x$CIbias$datastructure)
  {if (type=="estimate") {
    SSDvars <- c("SSD", "MMR", "BDI", "ERM", "FMA", "MoM", "BFM")
    CVvars <- c("CV", "CVsex")
    sdlogvars <- c("sdlog", "sdlogsex")
  }
  else if (type=="bias") {
    SSDvars <- c("MMR", "BDI", "ERM", "FMA", "MoM", "BFM")
    CVvars <- c("CV")
    sdlogvars <- c("sdlog")
  }}
  nSSDv <- sum(unique(x$estimates$methodUni) %in% SSDvars)
  nCVv <- sum(unique(x$estimates$methodUni) %in% CVvars)
  nsdlogv <- sum(unique(x$estimates$methodUni) %in% sdlogvars)
  nmethtotal <- nlevels(x$estimates$method)
  fntsize <- 11
  if (nmethtotal > 20) fntsize=9
  if (nmethtotal > 30) fntsize=8
  if (nmethtotal > 45) fntsize=7
  pltSSD <- NULL
  pltCV <- NULL
  pltsdlog <- NULL
  # generate pltSSD
  if (nSSDv > 0 & sum(is.na(x$estimates[x$estimates$methodUni %in% SSDvars,type])) < length(x$estimates[x$estimates$methodUni %in% SSDvars,type])) {
    {if (type=="estimate") {
	  pltSSD <- ggplot2::ggplot()
	  {if (attr(x$estimates, "estvalues")=="raw") pltSSD <- pltSSD + ggplot2::scale_x_log10() + ggplot2::geom_vline(xintercept=1) + ggplot2::xlab("estimate (axis scale is log10 transformed)")
	  else pltSSD <- pltSSD + ggplot2::geom_vline(xintercept=0)}
	  pltSSD <- pltSSD + ggplot2::geom_violin(data=x$estimates[x$estimates$methodUni %in% SSDvars,],
                          ggplot2::aes(y=method, x=estimate, fill=methodUni)) +
            ggplot2::stat_summary(data=x$estimates[x$estimates$methodUni %in% SSDvars,],
                          ggplot2::aes(y=method, x=estimate, fill=methodUni),fun=median, geom="point", size=2, color="black") +
            ggplot2::theme(legend.position = "none", axis.text.y=ggplot2::element_text(size=fntsize))
    }
	else if (type=="bias") pltSSD <- ggplot2::ggplot() + ggplot2::geom_vline(xintercept=0) +
            ggplot2::geom_violin(data=x$estimates[x$estimates$methodUni %in% SSDvars,],
                          ggplot2::aes(y=method, x=bias, fill=methodUni)) +
            ggplot2::stat_summary(data=x$estimates[x$estimates$methodUni %in% SSDvars,],
                          ggplot2::aes(y=method, x=bias, fill=methodUni),fun=median, geom="point", size=2, color="black") +
            ggplot2::theme(legend.position = "none", axis.text.y=ggplot2::element_text(size=fntsize)) +
			ggplot2::xlab("bias from sample SSD")}

    CIdat <- NULL
    if (type=="estimate" & !is.null(x$CI)) CIdat <- x$CI		  
    else if (type=="bias" & !is.null(x$CIbias)) CIdat <- x$CIbias			  
    if (!is.null(CIdat)) {		  
      # only do one or the other if upper or lower bound
	  segdatSSDlow <- NULL
	  segdatSSDhigh <- NULL
	  if (attr(x, "alternative")=="two.sided" | attr(x, "alternative")=="greater") {
        segdatSSDlow <- data.frame(x=CIdat$lower_lim[CIdat$methodUni %in% SSDvars],
                                 y=CIheight,
								 xend=CIdat$lower_lim[CIdat$methodUni %in% SSDvars],
								 yend=-CIheight)
        segdatSSDlow$y <- segdatSSDlow$y + as.integer(rownames(segdatSSDlow))
        segdatSSDlow$yend <- segdatSSDlow$yend + as.integer(rownames(segdatSSDlow))
	  }
	  if (attr(x, "alternative")=="two.sided" | attr(x, "alternative")=="less") {
        segdatSSDhigh <- data.frame(x=CIdat$upper_lim[CIdat$methodUni %in% SSDvars],
                                 y=CIheight,
								 xend=CIdat$upper_lim[CIdat$methodUni %in% SSDvars],
								 yend=-CIheight)
        segdatSSDhigh$y <- segdatSSDhigh$y + as.integer(rownames(segdatSSDhigh))
        segdatSSDhigh$yend <- segdatSSDhigh$yend + as.integer(rownames(segdatSSDhigh))
	  }
	  segdatSSD <- rbind(segdatSSDlow, segdatSSDhigh)
      pltSSD <- pltSSD + ggplot2::geom_segment(data=segdatSSD, ggplot2::aes(x=x, y=y, xend=xend, yend=yend),
                                       color="black")
      pltSSD <- pltSSD + ggplot2::geom_segment(data=segdatSSD, ggplot2::aes(x=x, y=y, xend=xend, yend=yend),
                                       color="red", linetype=3)
    }
  } # end pltSSD
  # generate pltCV
  if (nCVv > 0 & sum(is.na(x$estimates[x$estimates$methodUni %in% CVvars,type])) < length(x$estimates[x$estimates$methodUni %in% CVvars,type])) {
    {if (type=="estimate") pltCV <- ggplot2::ggplot() + ggplot2::geom_vline(xintercept=0) +
            ggplot2::geom_violin(data=x$estimates[x$estimates$methodUni %in% CVvars,],
                          ggplot2::aes(y=method, x=estimate, fill=methodUni)) +
            ggplot2::stat_summary(data=x$estimates[x$estimates$methodUni %in% CVvars,],
                          ggplot2::aes(y=method, x=estimate, fill=methodUni),fun=median, geom="point", size=2, color="black") +
            ggplot2::theme(legend.position = "none", axis.text.y=ggplot2::element_text(size=fntsize))
    else if (type=="bias") pltCV <- ggplot2::ggplot() + ggplot2::geom_vline(xintercept=0) +
            ggplot2::geom_violin(data=x$estimates[x$estimates$methodUni %in% CVvars,],
                          ggplot2::aes(y=method, x=bias, fill=methodUni)) +
            ggplot2::stat_summary(data=x$estimates[x$estimates$methodUni %in% CVvars,],
                          ggplot2::aes(y=method, x=bias, fill=methodUni),fun=median, geom="point", size=2, color="black") +
            ggplot2::theme(legend.position = "none", axis.text.y=ggplot2::element_text(size=fntsize)) +
			ggplot2::xlab("bias from sample CV")}
    CIdat <- NULL
    if (type=="estimate" & !is.null(x$CI)) CIdat <- x$CI		  
    else if (type=="bias" & !is.null(x$CIbias)) CIdat <- x$CIbias			  
    if (!is.null(CIdat)) {		  
      # only do one or the other if upper or lower bound
	  segdatCVlow <- NULL
	  segdatCVhigh <- NULL
	  if (attr(x, "alternative")=="two.sided" | attr(x, "alternative")=="greater") {
        segdatCVlow <- data.frame(x=CIdat$lower_lim[CIdat$methodUni %in% CVvars],
                                 y=CIheight,
								 xend=CIdat$lower_lim[CIdat$methodUni %in% CVvars],
								 yend=-CIheight)
        segdatCVlow$y <- segdatCVlow$y + as.integer(rownames(segdatCVlow))
        segdatCVlow$yend <- segdatCVlow$yend + as.integer(rownames(segdatCVlow))
	  }
	  if (attr(x, "alternative")=="two.sided" | attr(x, "alternative")=="less") {
        segdatCVhigh <- data.frame(x=CIdat$upper_lim[CIdat$methodUni %in% CVvars],
                                 y=CIheight,
								 xend=CIdat$upper_lim[CIdat$methodUni %in% CVvars],
								 yend=-CIheight)
        segdatCVhigh$y <- segdatCVhigh$y + as.integer(rownames(segdatCVhigh))
        segdatCVhigh$yend <- segdatCVhigh$yend + as.integer(rownames(segdatCVhigh))
	  }
	  segdatCV <- rbind(segdatCVlow, segdatCVhigh)
      pltCV <- pltCV + ggplot2::geom_segment(data=segdatCV, ggplot2::aes(x=x, y=y, xend=xend, yend=yend),
                                       color="black")
      pltCV <- pltCV + ggplot2::geom_segment(data=segdatCV, ggplot2::aes(x=x, y=y, xend=xend, yend=yend),
                                       color="red", linetype=3)
    }
  } # end pltCV
  # generate pltsdlog
  if (nsdlogv > 0 & sum(is.na(x$estimates[x$estimates$methodUni %in% sdlogvars,type])) < length(x$estimates[x$estimates$methodUni %in% sdlogvars,type])) {
    {if (type=="estimate") pltsdlog <- ggplot2::ggplot() + ggplot2::geom_vline(xintercept=0) +
            ggplot2::geom_violin(data=x$estimates[x$estimates$methodUni %in% sdlogvars,],
                          ggplot2::aes(y=method, x=estimate, fill=methodUni)) +
            ggplot2::stat_summary(data=x$estimates[x$estimates$methodUni %in% sdlogvars,],
                          ggplot2::aes(y=method, x=estimate, fill=methodUni),fun=median, geom="point", size=2, color="black") +
            ggplot2::theme(legend.position = "none", axis.text.y=ggplot2::element_text(size=fntsize))
    else if (type=="bias") pltsdlog <- ggplot2::ggplot() + ggplot2::geom_vline(xintercept=0) +
            ggplot2::geom_violin(data=x$estimates[x$estimates$methodUni %in% sdlogvars,],
                          ggplot2::aes(y=method, x=bias, fill=methodUni)) +
            ggplot2::stat_summary(data=x$estimates[x$estimates$methodUni %in% sdlogvars,],
                          ggplot2::aes(y=method, x=bias, fill=methodUni),fun=median, geom="point", size=2, color="black") +
            ggplot2::theme(legend.position = "none", axis.text.y=ggplot2::element_text(size=fntsize)) +
			ggplot2::xlab("bias from sample sdlog")}
    CIdat <- NULL
    if (type=="estimate" & !is.null(x$CI)) CIdat <- x$CI		  
    else if (type=="bias" & !is.null(x$CIbias)) CIdat <- x$CIbias			  
    if (!is.null(CIdat)) {		  
      # only do one or the other if upper or lower bound
	  segdatsdloglow <- NULL
	  segdatsdloghigh <- NULL
	  if (attr(x, "alternative")=="two.sided" | attr(x, "alternative")=="greater") {
        segdatsdloglow <- data.frame(x=CIdat$lower_lim[CIdat$methodUni %in% sdlogvars],
                                 y=CIheight,
								 xend=CIdat$lower_lim[CIdat$methodUni %in% sdlogvars],
								 yend=-CIheight)
        segdatsdloglow$y <- segdatsdloglow$y + as.integer(rownames(segdatsdloglow))
        segdatsdloglow$yend <- segdatsdloglow$yend + as.integer(rownames(segdatsdloglow))
	  }
	  if (attr(x, "alternative")=="two.sided" | attr(x, "alternative")=="less") {
        segdatsdloghigh <- data.frame(x=CIdat$upper_lim[CIdat$methodUni %in% sdlogvars],
                                 y=CIheight,
								 xend=CIdat$upper_lim[CIdat$methodUni %in% sdlogvars],
								 yend=-CIheight)
        segdatsdloghigh$y <- segdatsdloghigh$y + as.integer(rownames(segdatsdloghigh))
        segdatsdloghigh$yend <- segdatsdloghigh$yend + as.integer(rownames(segdatsdloghigh))
	  }
	  segdatsdlog <- rbind(segdatsdloglow, segdatsdloghigh)
      pltsdlog <- pltsdlog + ggplot2::geom_segment(data=segdatsdlog, ggplot2::aes(x=x, y=y, xend=xend, yend=yend),
                                       color="black")
      pltsdlog <- pltsdlog + ggplot2::geom_segment(data=segdatsdlog, ggplot2::aes(x=x, y=y, xend=xend, yend=yend),
                                       color="red", linetype=3)
    }
  } # end pltsdlog
  # build plotting text
  n <- nrow(x$sampleADS$comparative)
  nvar <- ncol(x$sampleADS$comparative)
  nResamp <- x$sampleADS$nResamp
  toptxt <- "Resampled"
  if (!is.null(attr(x, "resampling"))) {
    if (attr(x, "resampling")=="bootstrap") toptxt <- "Bootstrapped"
  }
  {if (type=="estimate") toptxt <- paste0(toptxt, " estimates")
  else if (type=="bias") toptxt <- paste0(toptxt, " estimate bias")}
  if (!is.null(attr(x, "conf.level"))) {
    toptxt <- paste0(toptxt, " with ")
	if (attr(x, "alternative")=="two.sided") toptxt <- paste0(toptxt, "two-sided ",
	     attr(x, "conf.level")*100, "% confidence intervals\n(", nvar, " variables, n=", n, ", n resampled datasets=", nResamp, ")")
	else if (attr(x, "alternative")=="less") toptxt <- paste0(toptxt, "one-sided ",
	     attr(x, "conf.level")*100, "% confidence intervals\n(upper bound, ", nvar, " variables, n=", n, ", n resampled datasets=", nResamp, ")")
	else if (attr(x, "alternative")=="greater") toptxt <- paste0(toptxt, "one-sided ", 
	     attr(x, "conf.level")*100, "% confidence intervals\n(lower bound, ", nvar, " variables, n=", n, ", n resampled datasets=", nResamp, ")")
  }
  {if (is.null(pltCV) & is.null(pltsdlog)) {
    return(gridExtra::grid.arrange(pltSSD, top=toptxt))
  }
  else if (is.null(pltSSD) & is.null(pltsdlog)) {
    return(gridExtra::grid.arrange(pltCV, top=toptxt))
  }
  else if (is.null(pltSSD) & is.null(pltCV)) {
    return(gridExtra::grid.arrange(pltsdlog, top=toptxt))
  }
  else if (is.null(pltSSD)) {
    return(gridExtra::grid.arrange(pltCV, pltsdlog,
	         heights=c(nCVv,nsdlogv), ncol = 1,
			 top=toptxt))
  }
  else if (is.null(pltCV)) {
    return(gridExtra::grid.arrange(pltSSD, pltsdlog,
	         heights=c(nSSDv,nsdlogv), ncol = 1,
			 top=toptxt))
  }
  else if (is.null(pltsdlog)) {
    return(gridExtra::grid.arrange(pltSSD, pltCV,
	         heights=c(nSSDv,nCVv), ncol = 1,
			 top=toptxt))
  }
  else {
      return(gridExtra::grid.arrange(pltSSD, pltCV, pltsdlog,
	         heights=c(nSSDv,nCVv,nsdlogv), ncol = 1,
			 top=toptxt))
  }}
}

#' @export
print.dimorphResampledMulti <- function(x, verbose=F) {
  methsUni <- levels(droplevels(x$estimates$methodUni))
  methsMulti <- levels(droplevels(x$estimates$methodMulti))
  centers <- levels(droplevels(x$estimates$center))
  centers[centers %in% "mean"] <- "arithmetic mean"
  centers[centers %in% "geomean"] <- "geometric mean"
  datastrucs <- levels(droplevels(x$estimates$datastructure))
  methcombos <- levels(droplevels(x$estimates$methodUni:x$estimates$methodMulti:x$estimates$center:x$estimates$datastructure))
  nmethtotal <- length(methcombos)
  txt <- "        dimorphResampledMulti Object\n\n"
  txt <- paste0(txt, "Comparative data set:\n")
  samplesize <- function(y) {
    nF <- 0
	nM <- 0
	nUnspecified <- 0
    if (is.null(y$compsex)) nUnspecified <- nrow(y$comparative)
	else {
	  nF <- sum(y$compsex==levels(y$compsex)[y$sex.female])
	  nM <- sum(y$compsex==levels(y$compsex)[(1:2)[-y$sex.female]])
	}
	out <- c(nFemale=nF, nMale=nM, nUnspecified=nUnspecified)
	return(out)
  }
  nsample <- samplesize(x$sampleADS)
  if (nsample["nFemale"]==0 & nsample["nMale"]==0) sextxt <- paste0(nsample["nUnspecified"], " (sex unspecified)")
  else {
    sextxt <- NULL
    if (nsample["nFemale"] > 0) sextxt <- paste0(sextxt, nsample["nFemale"], " female")
    if (nsample["nMale"] > 0) {
	  if (!is.null(sextxt)) sextxt <- paste0(sextxt, ", ")
	  sextxt <- paste0(sextxt, nsample["nMale"], " male")
	}
    if (nsample["nUnspecified"] > 0) {
	  if (!is.null(sextxt)) sextxt <- paste0(sextxt, ", ")
	  sextxt <- paste0(sextxt, nsample["nUnspecified"], " of unspecified sex")
	}
  }
  txt <- paste0(txt, "  number of specimens: ", sextxt, "\n")
  #txt <- paste0(txt, "  number of specimens: ", dim(x$sampleADS$comparative)[1], "\n")
  txt <- paste0(txt, "  number of variables: ", dim(x$sampleADS$comparative)[2], "\n")
  txt <- paste0(txt, "  variable names: ", paste0(colnames(x$sampleADS$comparative), collapse=", "), "\n")
  if ("TM" %in% methsMulti) txt <- paste0(txt, "  template variable: ", attr(x,"templatevar"), "\n")
  txt <- paste0(txt, "SSD estimate methods (univariate):\n  ", paste0(methsUni, collapse=", "), "\n")
  txt <- paste0(txt, "SSD estimate methods (multivariate):\n  ", paste0(methsMulti, collapse=", "), "\n")
  txt <- paste0(txt, "Centering algorithms:\n  ", paste0(centers, collapse=", "), "\n")
  txt <- paste0(txt, "Multivariate sampling with complete or missing data:\n  ", paste0(datastrucs, collapse=" and "), "\n")
  txt <- paste0(txt, "Number of unique combinations of univariate method, multivariate method,\n    centering algorithm, and complete or missing data structure: ", nmethtotal, "\n")
  if (verbose) txt <- paste0(txt, "Unique combinations:\n  ", paste0(methcombos, collapse="\n  "),"\n")
  txt <- paste0(txt, "Resampling data structure:\n")
  if (!is.null(x$sampleADS$exact)) {
    {if (x$sampleADS$exact) txt <- paste0(txt, "  type of resampling: exact")
    else if (!x$sampleADS$exact) txt <- paste0(txt, "  type of resampling: Monte Carlo")}
	if (!is.null(x$sampleADS$rebootstrapped)) if (x$sampleADS$rebootstrapped) txt <- paste0(txt, " for initial resampling, then Monte Carlo\n    for additional bootstrapping within subsamples")
	txt <- paste0(txt, "\n")
  }
  txt <- paste0(txt, "  number of resampled data sets: ", x$sampleADS$nResamp, "\n")
  strucClass <- x$sampleADS$strucClass
  if ("missing" %in% datastrucs) {
    structxt <- NULL
    {if (strucClass=="data.frame") structxt <- "\n    sampling individuals, then imposing missing data pattern"
    else if (strucClass=="vector") structxt <- "\n    sampling variable-specific number of individuals for each variable"}
    txt <- paste0(txt, "  missing data resampling structure: ", structxt, "\n")
  }
  {if (strucClass=="data.frame") {
    txt <- paste0(txt, "  number of individuals in each resampled data set: ", nrow(x$sampleADS$struc), "\n")
	txt <- paste0(txt, "  proportion of missing data in resampling structure: ", round(sum(is.na(x$sampleADS$struc)) / length(unlist(x$sampleADS$struc)), 3), "\n")
  }
  else if (strucClass=="vector") {
    txt <- paste0(txt, "  number of resampled individuals by variable in each resampled data set:\n",
	              paste0("    ", names(x$sampleADS$struc), ": ", x$sampleADS$struc, "\n", collapse=""))
  }}
  if (!is.null(attr(x, "resampling"))) txt <- paste0(txt, "  resampling procedure: ", attr(x, "resampling"), "\n")
  if (!is.null(x$sampleADS$replace)) {
    replacetxt <- NULL
    {if (x$sampleADS$replace) replacetxt <- "  subsamples sampled WITH replacement"
    else if (!x$sampleADS$replace) replacetxt <- "  subsamples sampled WITHOUT replacement"}
	if (!is.null(x$sampleADS$rebootstrapped)) if (x$sampleADS$rebootstrapped) replacetxt <- paste0(replacetxt, " for initial resampling, then sampled\n    WITH replacement within subsamples for additional bootstrapping")	
	txt <- paste0(txt, replacetxt, "\n")
  }
  if (!is.null(attr(x, "conf.level"))) {
    if (attr(x, "alternative")=="two.sided") txt <- paste0(txt, "  confidence intervals: two-sided, ", attr(x, "conf.level")*100, "% confidence level\n")
    if (attr(x, "alternative")=="less") txt <- paste0(txt, "  confidence intervals: one-sided (upper bound), ", attr(x, "conf.level")*100, "% confidence level\n")
    if (attr(x, "alternative")=="greater") txt <- paste0(txt, "  confidence intervals: one-sided (lower bound), ", attr(x, "conf.level")*100, "% confidence level\n")
  }
  {if (is.null(x$sampleADS$compsex)) sextxt <- "sex data absent"
  else sextxt <- "sex data present"}
  txt <- paste0(txt, "  other resampling parameters:\n    ", sextxt)  
  if (!is.null(attr(x$estimates, "estvalues"))) {
    if (attr(x$estimates, "estvalues")=="raw") logtxt <- "ratio variables (if present): ratio not logged"
	else if (attr(x$estimates, "estvalues")=="logged") logtxt <- "ratio variables (if present): natural log of ratio"
	txt <- paste0(txt, "\n    ", logtxt)
  }
  if (!is.null(attr(x, "matchvars"))) {
    if (attr(x, "matchvars")) mvtxt <- "matchvars = TRUE"
	else mvtxt <- "matchvars = FALSE"
	txt <- paste0(txt, "\n    ", mvtxt)
  }
  if (!is.null(attr(x, "na.rm"))) {
    if (attr(x, "na.rm")) narmtxt <- "na.rm = TRUE"
	else narmtxt <- "na.rm = FALSE"
	txt <- paste0(txt, "\n    ", narmtxt)
  }
  txt <- paste0(txt, "\n")
  cat(txt)
  # Add printing both the CI and CIbias data frames if present
  if (!is.null(x$CI)) {
    cat("\nConfidence intervals for estimates:\n")
    print(x$CI)
  }
  if (!is.null(x$CIbias)) {
    cat("\nConfidence intervals for bias of estimates from sample SSD, CVsex, or sdlogsex:\n")
    print(x$CIbias)
  }
}

#' @export
print.dimorphResampledUni <- function(x, verbose=F) {
  methsUni <- levels(droplevels(x$estimates$methodUni))
  centers <- levels(droplevels(x$estimates$center))
  centers[centers %in% "mean"] <- "arithmetic mean"
  centers[centers %in% "geomean"] <- "geometric mean"
  methcombos <- levels(droplevels(x$estimates$methodUni:x$estimates$center))
  nmethtotal <- length(methcombos)
  txt <- "        dimorphResampledUni Object\n\n"
  txt <- paste0(txt, "Comparative data set:\n")
  samplesize <- function(y) {
    nF <- 0
	nM <- 0
	nUnspecified <- 0
    if (is.null(y$compsex)) nUnspecified <- nrow(y$comparative)
	else {
	  nF <- sum(y$compsex==levels(y$compsex)[y$sex.female])
	  nM <- sum(y$compsex==levels(y$compsex)[(1:2)[-y$sex.female]])
	}
	out <- c(nFemale=nF, nMale=nM, nUnspecified=nUnspecified)
	return(out)
  }
  nsample <- samplesize(x$sampleADS)
  if (nsample["nFemale"]==0 & nsample["nMale"]==0) sextxt <- paste0(nsample["nUnspecified"], " (sex unspecified)")
  else {
    sextxt <- NULL
    if (nsample["nFemale"] > 0) sextxt <- paste0(sextxt, nsample["nFemale"], " female")
    if (nsample["nMale"] > 0) {
	  if (!is.null(sextxt)) sextxt <- paste0(sextxt, ", ")
	  sextxt <- paste0(sextxt, nsample["nMale"], " male")
	}
    if (nsample["nUnspecified"] > 0) {
	  if (!is.null(sextxt)) sextxt <- paste0(sextxt, ", ")
	  sextxt <- paste0(sextxt, nsample["nUnspecified"], " of unspecified sex")
	}
  }
  txt <- paste0(txt, "  number of specimens: ", sextxt, "\n")
  txt <- paste0(txt, "  number of variables: 1\n")
  txt <- paste0(txt, "  variable name: ", paste0(colnames(x$sampleADS$comparative), collapse=", "), "\n")
  txt <- paste0(txt, "SSD estimate methods (univariate):\n  ", paste0(methsUni, collapse=", "), "\n")
  txt <- paste0(txt, "Centering algorithms:\n  ", paste0(centers, collapse=", "), "\n")
  txt <- paste0(txt, "Number of unique combinations of univariate method and centering algorithm: ", nmethtotal, "\n")
  if (verbose) txt <- paste0(txt, "Unique combinations:\n  ", paste0(methcombos, collapse="\n  "),"\n")
  txt <- paste0(txt, "Resampling data structure:\n")
  if (!is.null(x$sampleADS$exact)) {
    {if (x$sampleADS$exact) txt <- paste0(txt, "  type of resampling: exact")
    else if (!x$sampleADS$exact) txt <- paste0(txt, "  type of resampling: Monte Carlo")}
	if (!is.null(x$sampleADS$rebootstrapped)) if (x$sampleADS$rebootstrapped) txt <- paste0(txt, " for initial resampling, then Monte Carlo\n    for additional bootstrapping within subsamples")
	txt <- paste0(txt, "\n")
  }
  txt <- paste0(txt, "  number of resampled data sets: ", x$sampleADS$nResamp, "\n")
  txt <- paste0(txt, "  number of individuals in each resampled data set: ", nrow(x$sampleADS$struc), "\n")
  if (!is.null(attr(x, "resampling"))) txt <- paste0(txt, "  resampling procedure: ", attr(x, "resampling"), "\n")
  if (!is.null(x$sampleADS$replace)) {
    replacetxt <- NULL
    {if (x$sampleADS$replace) replacetxt <- "  subsamples sampled WITH replacement"
    else if (!x$sampleADS$replace) replacetxt <- "  subsamples sampled WITHOUT replacement"}
	if (!is.null(x$sampleADS$rebootstrapped)) if (x$sampleADS$rebootstrapped) replacetxt <- paste0(replacetxt, " for initial resampling, then sampled\n    WITH replacement within subsamples for additional bootstrapping")	
	txt <- paste0(txt, replacetxt, "\n")
  }
  if (!is.null(attr(x, "conf.level"))) {
    if (attr(x, "alternative")=="two.sided") txt <- paste0(txt, "  confidence intervals: two-sided, ", attr(x, "conf.level")*100, "% confidence level\n")
    if (attr(x, "alternative")=="less") txt <- paste0(txt, "  confidence intervals: one-sided (upper bound), ", attr(x, "conf.level")*100, "% confidence level\n")
    if (attr(x, "alternative")=="greater") txt <- paste0(txt, "  confidence intervals: one-sided (lower bound), ", attr(x, "conf.level")*100, "% confidence level\n")
  }
  {if (is.null(x$sampleADS$compsex)) sextxt <- "sex data absent"
  else sextxt <- "sex data present"}
  txt <- paste0(txt, "  other resampling parameters:\n    ", sextxt)  
  if (!is.null(attr(x$estimates, "estvalues"))) {
    if (attr(x$estimates, "estvalues")=="raw") logtxt <- "ratio variables (if present): ratio not logged"
	else if (attr(x$estimates, "estvalues")=="logged") logtxt <- "ratio variables (if present): natural log of ratio"
	txt <- paste0(txt, "\n    ", logtxt)
  }
  if (!is.null(attr(x, "matchvars"))) {
    if (attr(x, "matchvars")) mvtxt <- "matchvars = TRUE"
	else mvtxt <- "matchvars = FALSE"
	txt <- paste0(txt, "\n    ", mvtxt)
  }
  if (!is.null(attr(x, "na.rm"))) {
    if (attr(x, "na.rm")) narmtxt <- "na.rm = TRUE"
	else narmtxt <- "na.rm = FALSE"
	txt <- paste0(txt, "\n    ", narmtxt)
  }
  txt <- paste0(txt, "\n")
  cat(txt)
  # Add printing both the CI and CIbias data frames if present
  if (!is.null(x$CI)) {
    cat("\nConfidence intervals for estimates:\n")
    print(x$CI)
  }
  if (!is.null(x$CIbias)) {
    cat("\nConfidence intervals for bias of estimates from sample SSD, CVsex, or sdlogsex:\n")
    print(x$CIbias)
  }
}


#' @export
print.dimorphAds <- function(x, verbose=F) {
  #multivar <- F
  #if (dim(x$comparative)[2] > 1) multivar <- T
  txt <- "        dimorphAds Object\n\n"
  txt <- paste0(txt, "Comparative data set:\n")
  samplesize <- function(y) {
    nF <- 0
	nM <- 0
	nUnspecified <- 0
    if (is.null(y$compsex)) nUnspecified <- nrow(y$comparative)
	else {
	  nF <- sum(y$compsex==levels(y$compsex)[y$sex.female])
	  nM <- sum(y$compsex==levels(y$compsex)[(1:2)[-y$sex.female]])
	}
	out <- c(nFemale=nF, nMale=nM, nUnspecified=nUnspecified)
	return(out)
  }
  nsample <- samplesize(x)
  if (nsample["nFemale"]==0 & nsample["nMale"]==0) sextxt <- paste0(nsample["nUnspecified"], " (sex unspecified)")
  else {
    sextxt <- NULL
    if (nsample["nFemale"] > 0) sextxt <- paste0(sextxt, nsample["nFemale"], " female")
    if (nsample["nMale"] > 0) {
	  if (!is.null(sextxt)) sextxt <- paste0(sextxt, ", ")
	  sextxt <- paste0(sextxt, nsample["nMale"], " male")
	}
    if (nsample["nUnspecified"] > 0) {
	  if (!is.null(sextxt)) sextxt <- paste0(sextxt, ", ")
	  sextxt <- paste0(sextxt, nsample["nUnspecified"], " of unspecified sex")
	}
  }
  txt <- paste0(txt, "  number of specimens: ", sextxt, "\n")
  txt <- paste0(txt, "  number of variables: ", dim(x$comparative)[2], "\n")
  txt <- paste0(txt, "  variable names: ", paste0(colnames(x$comparative), collapse=", "), "\n")
  txt <- paste0(txt, "Resampling data structure:\n")
  if (!is.null(x$exact)) {
    {if (x$exact) txt <- paste0(txt, "  type of resampling: exact")
    else if (!x$exact) txt <- paste0(txt, "  type of resampling: Monte Carlo")}
	if (!is.null(x$rebootstrapped)) if (x$rebootstrapped) txt <- paste0(txt, " for initial resampling, then Monte Carlo\n    for additional bootstrapping within subsamples")
	txt <- paste0(txt, "\n")
  }
  txt <- paste0(txt, "  number of resampled data sets: ", x$nResamp, "\n")
  strucClass <- x$strucClass
  structxt <- NULL
  {if (strucClass=="data.frame") structxt <- "\n    sampling individuals, then imposing missing data pattern"
  else if (strucClass=="vector") structxt <- "\n    sampling variable-specific number of individuals for each variable"}
  txt <- paste0(txt, "  missing data resampling structure: ", structxt, "\n")
  {if (strucClass=="data.frame") {
    txt <- paste0(txt, "  number of individuals in each resampled data set: ", nrow(x$struc), "\n")
	txt <- paste0(txt, "  proportion of missing data in resampling structure: ", round(sum(is.na(x$struc)) / length(unlist(x$struc)), 3), "\n")
  }
  else if (strucClass=="vector") {
    txt <- paste0(txt, "  number of resampled individuals by variable in each resampled data set:\n",
	              paste0("    ", names(x$struc), ": ", x$struc, "\n", collapse=""))
  }}
  if (!is.null(attr(x, "resampling"))) txt <- paste0(txt, "  resampling procedure: ", attr(x, "resampling"), "\n")
  if (!is.null(x$replace)) {
    replacetxt <- NULL
    {if (x$replace) replacetxt <- "  subsamples sampled WITH replacement"
    else if (!x$replace) replacetxt <- "  subsamples sampled WITHOUT replacement"}
	if (!is.null(x$rebootstrapped)) if (x$rebootstrapped) replacetxt <- paste0(replacetxt, " for initial resampling, then sampled\n    WITH replacement within subsamples for additional bootstrapping")	
	txt <- paste0(txt, replacetxt, "\n")
  }
  #{if (attr(x, "matchvars")) mvtxt <- "matchvars = TRUE"
  #else mvtxt <- "matchvars = FALSE"}
  #{if (attr(x, "na.rm")) narmtxt <- "na.rm = TRUE"
  #else narmtxt <- "na.rm = FALSE"}
  {if (is.null(x$compsex)) sextxt <- "sex data absent"
  else sextxt <- "sex data present"}
  txt <- paste0(txt, "  other resampling parameters:\n    ", sextxt)
  if (!is.null(x$matchvars)) {
    if (x$matchvars) mvtxt <- "matchvars = TRUE"
	else mvtxt <- "matchvars = FALSE"
	txt <- paste0(txt, "\n    ", mvtxt)
  }
  txt <- paste0(txt, "\n")
  cat(txt)
}


#' @export
print.SSDtest <- function(x, verbose=F, central="median") {
  central <- match.arg(arg=central, choices=c("median", "mean"))
  groups <- names(x$estimates)
  methcombos <- x$methcombos
  methsUni <- unique(methcombos$methodUni)
  methsMulti <- unique(methcombos$methodMulti)
  multi <- T
  if (is.null(methsMulti)) multi <- F
  centers <- unique(methcombos$center)
  centers[centers %in% "mean"] <- "arithmetic mean"
  centers[centers %in% "geomean"] <- "geometric mean"
  datastrucs <- unique(methcombos$datastructure)
  nmethtotal <- nrow(methcombos)
  {if (multi) methcombos <- methcombos[,c("methodUni", "methodMulti", "center", "datastructure", paste0(central, ".", groups))]
  else methcombos <- methcombos[,c("methodUni", "center", paste0(central, ".", groups))]}
  methcombos[,paste0(central, ".", groups)] <- round(methcombos[,paste0(central, ".", groups)], 4)
  colnames(methcombos) <- gsub(paste0(central, "."), "", colnames(methcombos))
  rownames(methcombos) <- NULL
  samplesize <- function(y) {
    nF <- 0
	nM <- 0
	nUnspecified <- 0
    if (is.null(y$sampleADS$compsex)) nUnspecified <- nrow(y$sampleADS$comparative)
	else {
	  nF <- sum(y$sampleADS$compsex==levels(y$sampleADS$compsex)[y$sampleADS$sex.female])
	  nM <- sum(y$sampleADS$compsex==levels(y$sampleADS$compsex)[(1:2)[-y$sampleADS$sex.female]])
	}
	out <- c(nFemale=nF, nMale=nM, nUnspecified=nUnspecified)
	return(out)
  }
  sampleDF <- data.frame(sample=groups, do.call(rbind, lapply(x$estimates, samplesize)))
  colnames(sampleDF) <- c("sample", "n female", "n male", "n unspecified")
  resamplingparams <- function(y) {
    #strucClass <- y$sampleADS$strucClass
	{if (y$sampleADS$exact) exact <- "exact"
	else exact <- "Monte Carlo"}
	rebootstrapped <- F
	if(!is.null(y$sampleADS$rebootstrapped)) rebootstrapped <- y$sampleADS$rebootstrapped
	if (rebootstrapped) exact <- paste0(exact, " for initial resampling, then Monte Carlo for additional bootstrapping within subsamples")
	nResamp <- y$sampleADS$nResamp
	replacement <- attr(y, "replace")
	{if (replacement) rpl <- "with replacement"
	else rpl <- "without replacement"}
	if (rebootstrapped) rpl <- paste0(rpl, " for initial resampling, then with replacement for additional bootstrapping within subsamples")
	out <- c("n resamples"=nResamp, "resampling type"=exact, "sampling"=rpl)
	return(out)
  }
  resampleDF <- data.frame(sample=groups, do.call(rbind, lapply(x$estimates, resamplingparams)))
  colnames(resampleDF) <- c("sample", "n resampled data sets", "resampling type", "sampling")
  txt <- "        SSDtest Object\n\n"
  txt <- paste0(txt, "Comparative data set:\n")
  cat(txt)
  print.data.frame(sampleDF, row.names = FALSE)
  #txt <- "\n"
  txt <- ""
  txt <- paste0(txt, "  number of variables: ", dim(x$estimates[[1]]$sampleADS$comparative)[2], "\n")
  txt <- paste0(txt, "  variable names: ", paste0(colnames(x$estimates[[1]]$sampleADS$comparative), collapse=", "), "\n")
  if ("TM" %in% x$methcombos$methodMulti) txt <- 
    paste0(txt, "  template variable: ", attr(x$estimates[[1]], "templatevar"), "\n")
  txt <- paste0(txt, "SSD estimate methods (univariate):\n  ", paste0(methsUni, collapse=", "), "\n")
  if (multi) txt <- paste0(txt, "SSD estimate methods (multivariate):\n  ", paste0(methsMulti, collapse=", "), "\n")
  txt <- paste0(txt, "Centering algorithms:\n  ", paste0(centers, collapse=", "), "\n")
  {if (multi) {
    txt <- paste0(txt, "Multivariate sampling with complete or missing data:\n  ", paste0(datastrucs, collapse=" and "), "\n")
    txt <- paste0(txt, "Number of unique combinations of univariate method, multivariate method,\n    centering algorithm, and complete or missing data structure: ", nmethtotal, "\n")
  }
  else {
    txt <- paste0(txt, "Number of unique combinations of univariate method and centering algorithm: ", nmethtotal, "\n")
  }}
  txt <- paste0(txt, "Resampling data structure:\n")
  #if (!is.null(x$estimates[[1]]$sampleADS$exact)) {
  #  {if (x$estimates[[1]]$sampleADS$exact) txt <- paste0(txt, "  type of resampling: exact")
  #  else if (!x$estimates[[1]]$sampleADS$exact) txt <- paste0(txt, "  type of resampling: Monte Carlo")}
  #	if (!is.null(x$estimates[[1]]$sampleADS$rebootstrapped)) if (x$estimates[[1]]$sampleADS$rebootstrapped) txt <- paste0(txt, " for initial resampling, then Monte Carlo\n    for additional bootstrapping within subsamples")
  #	txt <- paste0(txt, "\n")
  #}
  #txt <- paste0(txt, "  number of resampled data sets: ", x$estimates[[1]]$sampleADS$nResamp, "\n")
  cat(txt)
  print.data.frame(resampleDF, row.names = FALSE)
  strucClass <- x$estimates[[1]]$sampleADS$strucClass
  #txt <- "\n"
  txt <- ""
  if ("missing" %in% datastrucs) {
    structxt <- NULL
    {if (strucClass=="data.frame") structxt <- "\n    sampling individuals, then imposing missing data pattern"
    else if (strucClass=="vector") structxt <- "\n    sampling variable-specific number of individuals for each variable"}
    txt <- paste0(txt, "  missing data resampling structure: ", structxt, "\n")
  }
  #####
  # Need to add tag to SSDtest to indicate whether complete datasets were bootstrapped at full sample sizes
  fullsamplesboot <- F
  if (!is.null(attr(x$estimates[[1]], "fullsamplesboot"))) fullsamplesboot <- attr(x$estimates[[1]], "fullsamplesboot")
  #####
  {if (strucClass=="data.frame") {
    if (fullsamplesboot) txt <- paste0(txt, "  number of individuals in each resampled data set: total sample size within each group\n")
	else txt <- paste0(txt, "  number of individuals in each resampled data set: ", nrow(x$estimates[[1]]$sampleADS$struc), "\n")
	if (multi) txt <- paste0(txt, "  proportion of missing data in resampling structure: ", round(sum(is.na(x$estimates[[1]]$sampleADS$struc)) / length(unlist(x$estimates[[1]]$sampleADS$struc)), 3), "\n")
  }
  else if (strucClass=="vector") {
    txt <- paste0(txt, "  number of resampled individuals by variable in each resampled data set:\n",
	              paste0("    ", names(x$estimates[[1]]$sampleADS$struc), ": ", x$estimates[[1]]$sampleADS$struc, "\n", collapse=""))
  }}
  if (!is.null(attr(x$estimates[[1]], "resampling"))) txt <- paste0(txt, "  resampling procedure: ", attr(x, "resampling"), "\n")
#  if (!is.null(x$estimates[[1]]$sampleADS$replace)) {
#    replacetxt <- NULL
#    {if (x$estimates[[1]]$sampleADS$replace) replacetxt <- "  subsamples sampled WITH replacement"
#    else if (!x$estimates[[1]]$sampleADS$replace) replacetxt <- "  subsamples sampled WITHOUT replacement"}
#	if (!is.null(x$estimates[[1]]$sampleADS$rebootstrapped)) if (x$estimates[[1]]$sampleADS$rebootstrapped) replacetxt <- paste0(replacetxt, " for initial resampling, then sampled\n    WITH replacement within subsamples for additional bootstrapping")	
#	txt <- paste0(txt, replacetxt, "\n")
#  }
  #{if (is.null(x$sampleADS$compsex)) sextxt <- "sex data absent"
  #else sextxt <- "sex data present"}
  #txt <- paste0(txt, "  other resampling parameters:\n    ", sextxt)  
  txt <- paste0(txt, "  other resampling parameters:")  
  if (!is.null(attr(x$estimates[[1]]$estimates, "estvalues"))) {
    if (attr(x$estimates[[1]]$estimates, "estvalues")=="raw") logtxt <- "ratio variables (if present): ratio not logged"
	else if (attr(x$estimates[[1]]$estimates, "estvalues")=="logged") logtxt <- "ratio variables (if present): natural log of ratio"
	txt <- paste0(txt, "\n    ", logtxt)
  }
  if (!is.null(attr(x$estimates[[1]], "matchvars"))) {
    if (attr(x$estimates[[1]], "matchvars")) mvtxt <- "matchvars = TRUE"
	else mvtxt <- "matchvars = FALSE"
	txt <- paste0(txt, "\n    ", mvtxt)
  }
  if (!is.null(attr(x$estimates[[1]], "na.rm"))) {
    if (attr(x$estimates[[1]], "na.rm")) narmtxt <- "na.rm = TRUE"
	else narmtxt <- "na.rm = FALSE"
	txt <- paste0(txt, "\n    ", narmtxt)
  }
  txt <- paste0(txt, "\n")
  cat(txt)
  if (central=="mean") cnt <- "Mean"
  if (central=="median") cnt <- "Median"
  cat(paste0("\n", cnt, " resampled estimates for each combination of methods:\n"))
  print.data.frame(methcombos, row.names = T)
  if (!is.null(attr(x$estimates[[1]], "conf.level"))) {
    ## Add printing both the CI and CIbias data frames if present
    CI <- lapply(x$estimates, function(y) {out <- do.call(rbind,
	               apply(y$CI[,c("lower_lim", "upper_lim")], 1, function(z) paste0(round(z, 4), collapse=" to "), simplify=F))
                   #rownames(out) <- apply(y$CI[,!colnames(y$CI) %in% c("lower_lim", "upper_lim")], 1, function(z) paste0(z, collapse=":"))
				   rownames(out) <- apply(y$CI[,colnames(y$CI) %in% c("methodUni", "methodMulti", "center", "datastructure")], 1, function(z) paste0(z, collapse=":"))
                   return(out)
                 })
	methcomboCI <- methcombos
	rownames(methcomboCI) <- apply(methcomboCI[,!colnames(methcomboCI) %in% names(CI)], 1, function(z) paste0(z, collapse=":"))
	for (nm in names(CI)) {
	  methcomboCI[rownames(CI[[nm]]), nm] <- CI[[nm]]
	}
    if (attr(x$estimates[[1]], "alternative")=="two.sided") sidetxt <- "\nTwo-sided "
    if (attr(x$estimates[[1]], "alternative")=="less") sidetxt <- "\nOne-sided (upper bound) "
    if (attr(x$estimates[[1]], "alternative")=="greater") sidetxt <- "\nOne-sided (lower bound) " 
	txt <- paste0(sidetxt, attr(x$estimates[[1]], "conf.level")*100, "% confidence intervals for bootstrapped estimates:\n")
	cat(txt)
	print(methcomboCI)
    #if (!is.null(x$CIbias)) {
    #  cat("\nConfidence intervals for bias of estimates from sample SSD, CVsex, or sdlogsex:\n")
    #  print(x$CIbias)
    #}
  }
  if (!is.null(x$H0)) {
  cat(paste0("\n", cnt, " resampled estimates:\n"))
  print(H0)  
  }
  cat("\np-values (one-sided; null: first sample less or equally dimorphic as second sample):\n")
  pone <- round(x$pvalues$p.onesided,4)
  print(pone)
  cat("\np-values (two-sided):\n")
  ptwo <- round(x$pvalues$p.twosided,4)
  print(ptwo)
}

#' Plot \code{SSDtest} Object
#' 
#' Plot the object output by \code{\link[dimorph]{SSDtest}}.  \strong{N.B.:} warning messages stating 
#'   "\code{Removed} \emph{n} \code{rows containing  missing values (`geom_bar()`).}" are generated 
#'   for any histogram that contains bars of height zero, which will be true for most histograms.  This is 
#'   due to a known issue with how \code{\link[ggplot2]{ggplot}} handles histograms (see 
#'   \href{https://github.com/tidyverse/ggplot2/issues/3265}{https://github.com/tidyverse/ggplot2/issues/3265}), 
#'   and these warning messages can safely be ignored.  
#' @param x An \code{SSDtest} object.
#' @param est An integer specifiying the specific combination of methods (e.g., univariate 
#'   method, centering algorithm) for the resampled estimates to be plotted.  This integer corresponds to 
#'   the row number in the table providing mean or median estimates by method combination when printing
#'   the \code{SSDtest} object.
#' @param type A character string specifying the type of plot to produce. \code{"est"} generates a plot
#'   of histograms of estimates for all of the samples in \code{x}. \code{"diff"} gereates a plot of one 
#'   or more histograms of all possible pairwise differences in estimates between pairs of samples, with 
#'   the argument \code{diffs} specifying which pairs to plot.  Defaults to \code{"est"}.  Also, because 
#'   the plotting function calculates histograms
#'   for \code{type="diff"} by first calculating the difference between all possible pairs of one estimate taken 
#'   from the two distributions of resampled values to be compared, plots using \code{type="diff"} can take a 
#'   long time to be generated.  The output of \code{plot.SSDtest()} can be saved and plotted later without the need 
#'   for recalculation.
#' @param diffs Specifies which set of samples to plot pairwise differences for if \code{type} is \code{"diff"}. If set
#'   to \code{"all"} it generates a single plot with all possible pairs.  If given a vector of two or more integers,
#'   it plots all of the differences for that set of samples (e.g., setting \code{diffs} to \code{c(1,3,4)} will plot 
#'   three histograms: one for the differences in estimates between the first and third sample in \code{x}, one for the
#'   differences between the first and fourth sample, and one for the differences between the third and fourth sample).
#'   Any of the samples identified in \code{diffs} that do not include data for the estimate specified in \code{est} 
#'   will not be included in pairwise comparisons (e.g., if \code{"SSD"} was included as a univariate method on a set 
#'   of samples for which some had sex data and others didn't, the samples without sex data would be missing estimates 
#'   for that method and thus would be excluded from pairwise comparisons of \code{"SSD"} estimates.) An error will
#'   be generated in cases where fewer then two samples with estimates are specified.  Defaults to "all".
#' @param nbins An integer specifiying the number of bins to sort data into for histograms.
#' @param plottitle A character string specifying a title to print at the top of the plot.  If \code{NULL}, no 
#'   title is printed.  Defaults to \code{NULL}.
#' @param groupcols A vector of hexidecimal colors of length equal to the number of samples in \code{x} for use
#'   in plots when \code{type} is \code{"est"}.  Defaults to a palette generated by \code{\link[pals]{parula}}.
#' @param leg A logical scalar indicating whether a legend should be provided.  Defaults to \code{TRUE}. 
#' @param legpos A character string specifying where to plot the legend.  Acceptable values are \code{"left"}, 
#'   \code{"top"}, \code{"right"}, and \code{"bottom"}.  Defaults to \code{"right"}.
#' @param legsize A numeric scalar specifying a multiplier for the legend font size relative to the default font size. 
#'   Defaults to \code{1}.
#' @param legtitle A character string specifying a title to print at the top of the legend.  If \code{NULL}, no 
#'   title is printed.  Defaults to \code{NULL}.
#' @param titlesize A numeric scalar specifying a multiplier for the title font size relative to the default font size. 
#'   Defaults to \code{2}.
#' @param invert An integer vector specifying the group number of histograms to invert and project below the zero 
#'   count line, where group number is determined by the order of datasets provided first in \code{fossil}, then in 
#'   \code{comp}.  If \code{NULL}, no histograms are inverted.  Defaults to \code{NULL}.
#' @param xlim A numeric vector with two elements specifying the minimum and maximum values of the x-axis.  If 
#'   \code{xlim} is \code{NULL}, default values generated by \code{\link[ggplot2]{ggplot}} are used.
#' @param ylim A numeric vector with two elements specifying the minimum and maximum values of the y-axis.  If 
#'   \code{ylim} is \code{NULL}, default values generated by \code{\link[ggplot2]{ggplot}} are used.
#' @param diffcol A color used to fill in histogram bars for differences in parameter values between two samples. 
#'   defaults to \code{"#CCB03D"}.
#' @examples
#' SSDvars <- c("FHSI", "TPML", "TPMAP", "TPLAP", "HHMaj", "HHMin", "RHMaj", "RHMin", "RDAP", "RDML")
#' test_faux_multi1 <- SSDtest(fossil=list("Fauxil sp. 1"=fauxil[fauxil$Species=="Fauxil sp. 1", SSDvars]),
#'                             comp=list("G. gorilla"=apelimbart[apelimbart$Species=="Gorilla gorilla", SSDvars],
#'                                       "H. sapiens"=apelimbart[apelimbart$Species=="Homo sapiens", SSDvars],
#'                                       "P. troglodytes"=apelimbart[apelimbart$Species=="Pan troglodytes", SSDvars],
#'                                       "H. lar"=apelimbart[apelimbart$Species=="Hylobates lar", SSDvars]),
#'                             fossilsex=NULL,
#'                             compsex=list("G. gorilla"=apelimbart[apelimbart$Species=="Gorilla gorilla", "Sex"],
#'                                          "H. sapiens"=apelimbart[apelimbart$Species=="Homo sapiens", "Sex"],
#'                                          "P. troglodytes"=apelimbart[apelimbart$Species=="Pan troglodytes", "Sex"],
#'                                          "H. lar"=apelimbart[apelimbart$Species=="Hylobates lar", "Sex"]),
#'                             methsUni=c("SSD", "MMR", "BDI"),
#'                             methsMulti=c("GMM"),
#'                             datastruc="both",
#'                             nResamp=100,
#'                             templatevar="FHSI")
#' test_faux_multi1
#' plot(test_faux_multi1)
#' speciescolors <- c("Fauxil sp. 1"="#352A87", "G. gorilla"="#EABA4B",
#'                    "H. sapiens"="#09A9C0", "P. troglodytes"="#78BE7C",
#'                    "H. lar"="#0D77DA")
#' plot(test_faux_multi1, groupcols=speciescolors) # change colors of sample distributions
#' plot(test_faux_multi1, est=4, groupcols=speciescolors) # plot estimates for the fourth methcombo (GMM and MMR)
#' plot(test_faux_multi1, est=4, groupcols=speciescolors, invert=c(2,4)) # invert second and fourth sample distributions
#' plot(test_faux_multi1, est=4, type="diff") # plot distributions of differences among sample estimates for methocombo 4
#' plot(test_faux_multi1, est=4, type="diff", diffs=c(2,3)) # plots differences between second and third samples
#' @export
plot.SSDtest <- function(x, est=1, type="est", diffs=NULL, nbins=100, plottitle=NULL,
                         groupcols=NULL, leg=T, legpos=NULL, legsize=1, legtitle="", titlesize=2,
						 invert=NULL, xlim=NULL, ylim=NULL,
						 diffcol="#CCB03D") {
  type <- match.arg(type, choices=c("est", "diff"))
  ndat <- length(x$estimates)
  groups <- names(x$estimates)
  # 'est' is row number for methcombo providing the specific estimate to plot
  methcombos <- x$methcombos
  if (!est %in% 1:nrow(methcombos)) stop("'est' must be an integer specifying the row number of the method combination in 'x' to plot.")
  methsUni <- unique(methcombos$methodUni)
  methsMulti <- unique(methcombos$methodMulti)
  multi <- T
  if (is.null(methsMulti)) multi <- F
  centers <- unique(methcombos$center)
  centers[centers %in% "mean"] <- "arithmetic mean"
  centers[centers %in% "geomean"] <- "geometric mean"
  datastrucs <- unique(methcombos$datastructure)
  {if (type=="est") {
    if (is.null(groupcols)) {
      groupcols <- pals::parula(ndat)
	  names(groupcols) <- groups
    }
    {if (multi) {
	  dat <- lapply(x$estimates, function(y) y$estimates$estimate[y$estimates$methodUni==methcombos$methodUni[est] &
	                                                              y$estimates$methodMulti==methcombos$methodMulti[est] &
																  y$estimates$center==methcombos$center[est] &
																  y$estimates$datastructure==methcombos$datastructure[est]])
	}
	else {
	  dat <-  lapply(x$estimates, function(y) y$estimates$estimate[y$estimates$methodUni==methcombos$methodUni[est] &
																   y$estimates$center==methcombos$center[est]])
 	}}
	for (i in 1:ndat) dat[[i]] <- dat[[i]][!is.na(dat[[i]])]
	datrange <- range(dat, na.rm=T)
	if (!is.null(xlim)) {
	  datrange <- range(datrange, xlim)
	  if (!identical(datrange, xlim)) warning("The data range to plot extends beyond the user-specified 'xlim', so 'xlim' will be adjusted.")
	}
	xlim <- datrange
	# extend plot 1% to either side of data range
	xlim[1] <- xlim[1]-diff(datrange)*.01
	xlim[2] <- xlim[2]+diff(datrange)*.01
	bw <- round(diff(datrange)/100,3)
	if (bw < 0.001) bw <- 0.001
	# Combine into single data frame.  If any sample has only a single data point then hold it out separate
	#   to be plotted as a line. This occurs when only one sample has missing data, i.e., a fossil sample, and
	#   a point estimate is generated for that sample.
	plotdat <- NULL
	singletonest <- NULL
	singletongroup <- NULL
	for (i in 1:ndat) {
	  if (length(dat[[i]])< 2) {
	    if (length(dat[[i]])==0) next
	    #if (is.na(dat[[i]][1])) next
		singletonest <- dat[[i]][1]
		singletongroup <- groups[i]
		next
	  }
	  plotdat <- rbind(plotdat, data.frame(group=rep(groups[i], length(dat[[i]])), estimate=dat[[i]]))
	}
	plotdat$group <- factor(plotdat$group, levels=groups[groups %in% unique(plotdat$group)])
    datatxt <- NULL
	if (multi) {
	  if (methcombos$datastructure[est]=="complete") datatxt <- "(no missing data)"
      else if (methcombos$datastructure[est]=="missing") datatxt <- "(missing data)"
	}
	# Generate plot
	plt <- ggplot2::ggplot(plotdat, ggplot2::aes(x=estimate, fill=group)) + 
	         ggplot2::scale_fill_manual(values=groupcols[levels(plotdat$group)]) + 
			 ggplot2::guides(fill=ggplot2::guide_legend(title=legtitle)) +
			 ggplot2::scale_x_continuous(limits=xlim)
	if (!is.null(ylim)) plt <- plt + ggplot2::scale_y_continuous(limits=ylim)
	{if (multi) {
	  if (attr(x$estimates[[1]]$estimates, "estvalues")=="logged")
	    {if (methcombos$methodUni[est] %in% c("CV", "CVsex", "sdlog", "sdlogsex")) plt <- plt + ggplot2::labs(x=bquote(.(methcombos$methodUni[est])[.(methcombos$methodMulti[est])]~.(datatxt)))
	    else plt <- plt + ggplot2::labs(x=bquote(ln~.(methcombos$methodUni[est])[.(methcombos$methodMulti[est])]~.(datatxt)))}
	  else
	    plt <- plt + ggplot2::labs(x=bquote(.(methcombos$methodUni[est])[.(methcombos$methodMulti[est])]~.(datatxt)))
	}
	else { # univariate
	  if (attr(x$estimates[[1]]$estimates, "estvalues")=="logged")
	    {if (methcombos$methodUni[est] %in% c("CV", "CVsex", "sdlog", "sdlogsex")) plt <- plt + ggplot2::labs(x=bquote(.(methcombos$methodUni[est])))
	    else plt <- plt + ggplot2::labs(x=bquote(ln~.(methcombos$methodUni[est])))}
	  else
	    plt <- plt + ggplot2::labs(x=bquote(.(methcombos$methodUni[est])))
	}}
	# plot horizontal black line at zero
	plt <- plt + ggplot2::geom_hline(yintercept=0, col="black")
	{if (is.null(invert)) {
	  plt <- plt +
			 ggplot2::geom_histogram(ggplot2::aes(y=ggplot2::after_stat(density)), binwidth=bw, alpha=1, position="identity") +
			 ggplot2::geom_histogram(ggplot2::aes(y=ggplot2::after_stat(density)), binwidth=bw, alpha=.4, position="identity")
	}
	else {
	  datdown <- plotdat[plotdat$group %in% groups[invert],]
	  datup <- plotdat[plotdat$group %in% groups[-invert],]
	  plt <- plt +
			 ggplot2::geom_histogram(data=datup, ggplot2::aes(y=ggplot2::after_stat(density)), binwidth=bw, alpha=1, position="identity") +
			 ggplot2::geom_histogram(data=datup, ggplot2::aes(y=ggplot2::after_stat(density)), binwidth=bw, alpha=.4, position="identity") +
			 ggplot2::geom_histogram(data=datdown, ggplot2::aes(y=ggplot2::after_stat(density)*-1), binwidth=bw, alpha=1, position="identity") +
			 ggplot2::geom_histogram(data=datdown, ggplot2::aes(y=ggplot2::after_stat(density)*-1), binwidth=bw, alpha=.4, position="identity")
	}}
	if (!is.null(singletonest)) {
	  #plt <- plt + ggplot2::geom_vline(xintercept=singletonest, linetype="dashed", color="black", linewidth=1
	  plt <- plt + ggplot2::geom_vline(ggplot2::aes(xintercept=singletonest, color=singletongroup),
	                                   linetype="dashed", linewidth=1) +
				   ggplot2::scale_color_manual(name="", values = "black")
	}
	# add title
	if (!is.null(plottitle)) plt <- plt + ggplot2::ggtitle(plottitle)
    # plot legend
	{if (leg) {
	  if (is.null(legpos)) legpos <- "right"
	  plt <- plt + ggplot2::theme(legend.position=legpos,
	                        legend.text=ggplot2::element_text(face="italic", size=ggplot2::rel(legsize)),
							plot.title=ggplot2::element_text(face="bold", size=ggplot2::rel(titlesize)),
							legend.title=ggplot2::element_text(face="bold", size=ggplot2::rel(legsize)))
	}
	else plt <- plt + ggplot2::theme(legend.position="none",
	                                 plot.title=ggplot2::element_text(face="bold", size=ggplot2::rel(titlesize)))}
	# remove y-axis ticks and labels
	#plt <- plt + ggplot2::theme(axis.text.y=ggplot2::element_blank(), axis.ticks.y=ggplot2::element_blank()) 
  } # end 'est'
  ########################################################################
  else if (type=="diff") {
    # starts the same as 'est'
	{if (multi) {
	  dat <- lapply(x$estimates, function(y) y$estimates$estimate[y$estimates$methodUni==methcombos$methodUni[est] &
	                                                              y$estimates$methodMulti==methcombos$methodMulti[est] &
																  y$estimates$center==methcombos$center[est] &
																  y$estimates$datastructure==methcombos$datastructure[est]])
	}
	else {
	  dat <-  lapply(x$estimates, function(y) y$estimates$estimate[y$estimates$methodUni==methcombos$methodUni[est] &
																   y$estimates$center==methcombos$center[est]])
 	}}
	for (i in 1:ndat) dat[[i]] <- dat[[i]][!is.na(dat[[i]])]
    # remove components of dat that are not called for by user
	if (is.null(diffs)) diffs <- "all"
	differrortxt <- "'diffs' must either be 'all' or a vector of integers specifying the samples to be included."
	{if (class(diffs)=="numeric"|class(diffs)=="integer") {
	  if(prod(diffs==as.integer(diffs))==0) stop(differrortxt)
      diffsnew <- diffs[diffs %in% 1:ndat]
	  if (length(diffsnew)==0) stop(differrortxt)
	  if (!identical(diffs, diffsnew)) warning("Some integer values in 'diffs' do not correspond to samples in 'x' and have been removed.")
	  diffs <- diffsnew
	  rm(diffsnew)
	  for (i in 1:ndat) {
	    if (!(groups[i] %in% groups[diffs])) dat[[groups[i]]] <- NULL
	  }
	}
	else if(!identical(diffs, "all")) stop(differrortxt)}
	ndatnew <- length(dat)
	groupsnew <- names(dat)
	# now needs to calculate sets of differences - pull code from SSDtest
	pair <- data.frame(sample.greater=rep(NA, choose(ndatnew,2)), sample.lesser=rep(NA, choose(ndatnew,2)))
	counter <- 0
    for (i in 1:(ndatnew-1)) {
	  for (j in (i+1):ndatnew) {
	    counter <- counter+1
		pair$sample.greater[counter] <- groupsnew[i]
		pair$sample.lesser[counter] <- groupsnew[j]
	  }
	}
	xDiffs <- apply(pair, 1, function(y) rep(dat[[as.character(y["sample.greater"])]],
	                                         each=length(dat[[as.character(y["sample.lesser"])]])) -
	                                     rep(dat[[as.character(y["sample.lesser"])]],
										     length(dat[[as.character(y["sample.greater"])]])),
					simplify=F)
	names(xDiffs) <- paste0(pair$sample.greater, " - ", pair$sample.lesser)
	dat <- xDiffs
	ncomparisons <- length(dat)
	comparisons <- names(dat)
	for (i in 1:ncomparisons) if (length(dat[[comparisons[i]]])==0) dat[[comparisons[i]]] <- NULL
	if (length(dat)==0) stop("The samples specified by 'diffs' do not include at least two samples with data for this value of 'est'.")
	ncomparisons <- length(dat)
	comparisons <- names(dat)
	datrange <- range(dat, na.rm=T)
	if (!is.null(xlim)) {
	  datrange <- range(datrange, xlim)
	  if (!identical(datrange, xlim)) warning("The data range to plot extends beyond the user-specified 'xlim', so 'xlim' will be adjusted.")
	}
	xlim <- datrange
	# extend plot 1% to either side of data range
	xlim[1] <- xlim[1]-diff(datrange)*.01
	xlim[2] <- xlim[2]+diff(datrange)*.01
	bw <- round(diff(datrange)/nbins,3)
	if (bw < 0.001) bw <- 0.001
	rm(xDiffs)
	datatxt <- NULL
	if (multi) {
	  if (methcombos$datastructure[est]=="complete") datatxt <- "(no missing data)"
      else if (methcombos$datastructure[est]=="missing") datatxt <- "(missing data)"
	}
	# Generate plot
	pltlist <- list(NULL)
	for (i in 1:ncomparisons) {
	  pvalone <- round(x$pvalues$p.onesided[,est][names(dat)[i]], 3)
	  if (pvalone==0) pvalone <- "< 0.001"
	  else pvalone <- paste0("= ", sprintf(pvalone, fmt = '%#.3f'))
	  pvaltwo <- round(x$pvalues$p.twosided[,est][names(dat)[i]], 3)
	  if (pvaltwo==0) pvaltwo <- "< 0.001"
	  else pvaltwo <- paste0("= ", sprintf(pvaltwo, fmt = '%#.3f'))
	  # generate ggplot in pltlist element 'i'
	  if (is.null(plottitle))
	  pltlist[[i]] <- ggplot2::ggplot(data.frame(difference=dat[[i]]), ggplot2::aes(x=difference)) + 
						#ggplot2::xlim(datrange) +
						#ggplot2::xlim(xlim) +
						ggplot2::geom_histogram(ggplot2::aes(y=ggplot2::after_stat(density)),
						                        binwidth=bw, alpha=1, position="identity", fill=diffcol)
      {if (is.null(plottitle)) pltlist[[i]] <- pltlist[[i]] + ggplot2::ggtitle(names(dat)[i],
						                 subtitle=bquote("one-sided" ~ italic(p)~.(pvalone)~"|"~"two-sided" ~ italic(p)~.(pvaltwo)))
	  else pltlist[[i]] <- pltlist[[i]] + ggplot2::ggtitle(plottitle)}
      pltlist[[i]] <- pltlist[[i]] + ggplot2::scale_x_continuous(limits=xlim)
	  if (!is.null(ylim)) pltlist[[i]] <- pltlist[[i]] + ggplot2::scale_y_continuous(limits=ylim)
	  # plot horizontal black line at zero
	  pltlist[[i]] <- pltlist[[i]] + ggplot2::geom_hline(yintercept=0, col="black")
	  {if (multi) {
	    if (attr(x$estimates[[1]]$estimates, "estvalues")=="logged")
	      {if (methcombos$methodUni[est] %in% c("CV", "CVsex", "sdlog", "sdlogsex")) pltlist[[i]] <- pltlist[[i]] + ggplot2::labs(x=bquote(Delta ~ .(methcombos$methodUni[est])[.(methcombos$methodMulti[est])]~.(datatxt)))
	      else pltlist[[i]] <- pltlist[[i]] + ggplot2::labs(x=bquote(Delta ~ ln~.(methcombos$methodUni[est])[.(methcombos$methodMulti[est])]~.(datatxt)))}
	    else
	      pltlist[[i]] <- pltlist[[i]] + ggplot2::labs(x=bquote(Delta ~ .(methcombos$methodUni[est])[.(methcombos$methodMulti[est])]~.(datatxt)))
	  }
	  else { # univariate
	    if (attr(x$estimates[[1]]$estimates, "estvalues")=="logged") 
	      {if (methcombos$methodUni[est] %in% c("CV", "CVsex", "sdlog", "sdlogsex")) pltlist[[i]] <- pltlist[[i]] + ggplot2::labs(x=bquote(Delta ~ .(methcombos$methodUni[est])))
	      else pltlist[[i]] <- pltlist[[i]] + ggplot2::labs(x=bquote(Delta ~ ln~.(methcombos$methodUni[est])))}
	    else
	      pltlist[[i]] <- pltlist[[i]] + ggplot2::labs(x=bquote(Delta ~ .(methcombos$methodUni[est])))
	  }}
      pltlist[[i]] <- pltlist[[i]] + ggplot2::geom_vline(xintercept=0, linetype="solid", color="black", linewidth=1)
	  # remove y-axis ticks and labels
	  #pltlist[[i]] <- pltlist[[i]] + ggplot2::theme(axis.text.y=ggplot2::element_blank(), axis.ticks.y=ggplot2::element_blank())   
	} # end for loop 'ncomparisons'
	# convert pltlist to plt
	{if (ncomparisons==1) {
	  plt <- pltlist[[1]]
	  rm(pltlist)
	}
	else {
	  ngridcol <- ceiling(sqrt(length(dat)))
	  ngridrow <- ceiling(length(dat)/ngridcol)
	  plt <- gridExtra::marrangeGrob(grobs=pltlist, nrow=ngridrow, ncol=ngridcol,
	           top="")
	}}
  } # end 'diff'
  } # end if else
  return(plt)
}
