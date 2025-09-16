# ============================== CODE METADATA =============================== #
# AUTHOR: Fernando Colchero
# DATE CREATED: 2025-10x-11
# DESCRIPTION: Additional functions to run multiple BaFTA analyses.
# NOTES: 
# ================================ CODE START ================================ #
# Multiple BaFTA analyses:
multibafta <- function(object, dataType = "aggregated", models, 
                       minAge = NA, gestTime = NA, 
                       niter = 55000, burnin = 5001, thinning = 20, 
                       nsim = 4, ncpus = 4) {
  allModels <- c("quadratic", "PeristeraKostaki", "ColcheroMuller", 
                 "Hadwiger", "gamma", "beta", "skewNormal", 'gammaMixture',
                 "HadwigerMixture", "skewSymmetric", "skewLogistic")
  modList <- paste(paste("'", allModels, "'", sep = ""), collapse = ", ")
  
  if (any(!models %in% allModels)) {
    stop(sprintf("Wrong model specification. Available models are: %s.\n", 
                 modList))
  }
  
  # Number of models:
  nmods <- length(models)
  
  # list of results:
  inflist <- list()
  gofdf <- data.frame(model = models, DIC = rep(NA, nmods), 
                      PredLoss = rep(NA, nmods), stringsAsFactors = FALSE)
  
  # run analyses:
  for (imod in 1:nmods) {
    # Print to the console the current model:
    updt <- sprintf("Run number %s, model: %s", imod, models[imod])
    cat(paste("\n", paste(rep("-", nchar(updt)), collapse = ""), sep = ""))
    cat(sprintf("\n%s\n", updt))
    cat(paste(paste(rep("-", nchar(updt)), collapse = ""), "\n", sep = ""))
    
    # run model:
    out <- bafta(object = object, dataType = dataType, model = models[imod], 
                 minAge = minAge, gestTime = gestTime, 
                 niter = niter, burnin = burnin, thinning = thinning, 
                 nsim = nsim, ncpus = ncpus)
    inflist[[models[imod]]] <- out
    if (all(out$coefficients[, "Rhat"] < 1.05)) {
      gofdf$DIC[imod] <- out$DIC["DIC"]
      gofdf$PredLoss[imod] <- out$PredLoss[1, "Deviance"]
    }
  }
  
  # Sort by lowest predictive loss:
  idsort <- list(PredLoss = sort.int(gofdf$PredLoss, index.return = TRUE, 
                                     na.last = TRUE)$ix,
                 DIC = sort.int(gofdf$DIC, index.return = TRUE, 
                                na.last = TRUE)$ix)
  
  # Create output list:
  outList <- list(runs = inflist, gof = gofdf, idsort = idsort)
  class(outList) <- "multibafta"
  return(outList)
}

# Plot for multi-BaFTA:
plot.multibafta <- function(x, sortBy = "DIC") {
  # All models in BaFTA:
  allModels <- c("quadratic", "PeristeraKostaki", "ColcheroMuller", 
                 "Hadwiger", "gamma", "beta", "skewNormal", 'gammaMixture',
                 "HadwigerMixture", "skewSymmetric", "skewLogistic")
  
  # Labels for all models:
  modelLabs <- c("Quadratic", "Peristera-Kostaki", "Colchero-Muller",
                 "Hadwiger", "Gamma", "Beta", "Skew-normal", "Gamma mixture",
                 "Hadwiger mixture", "Skew-symmetric", "Skew-logistic")
  names(modelLabs) <- allModels
  
  # number of models:
  nmods <- nrow(x$gof)
  
  # models ran:
  models <- x$gof$model
  
  # ran model names:
  modnames <- modelLabs[models]
  
  # Layout matrix:
  laymat <- rbind(cbind(rep(2, nmods * 2), seq(1, nmods * 4, 2) + 5, 
                        rep(3, nmods * 2), seq(2, nmods * 4, 2) + 5), 
                  c(0, 4, 0, 5), c(0, rep(1, 3)))
  
  # Heights and widths:
  heights <- c(rep(c(0.15, 0.5), nmods), 0.1, 0.2)
  widths <- rep(c(0.15, 1.5), 2)
  whr <- sum(widths) / sum(heights)
  
  # X and Y limits:
  xlim <- c(0, max(x$runs[[1]]$aggrData$Age) + 1)
  ylim <- list(fert = c(0, max(sapply(1:nmods, function(imod) {
    mf <- max(x$runs[[imod]]$fert[, "Upper"], na.rm = TRUE)
    if (mf == Inf) mf <- 0.1
    return(mf)
  }))),
  pred = c(0, max(sapply(1:nmods, function(imod) {
    mp <- max(x$runs[[imod]]$pred[, "Upper"], na.rm = TRUE)
    if (mp == Inf) mp <- 1
    return(mp)
  }))))
  
  # Margins:
  mar <- c(0.5, 2, 0.5, 0.5)
  
  
  # Produce plot:
  layout(mat = laymat, widths = widths, heights = heights)
  
  # X-axis label:
  par(mar = mar * c(0, 1, 0, 1))
  plot(c(0, 1), c(0, 1), col = NA, xlab = "", ylab = "", axes = FALSE)
  text(0.5, 0.95, "Age", cex = 2, xpd = NA)
  lines(c(0.1, 0.25), c(0.3, 0.3), type = 'b', lwd = 2)
  text(0.27, 0.3, "Observed", adj = 0, cex = 1.25)
  polygon(c(0.55, 0.7)[c(1, 2, 2, 1)], c(0.1, 0.5)[c(1, 1, 2, 2)], 
          col = 'orange', border = NA)
  lines(c(0.55, 0.7), c(0.3, 0.3), col = 'red', lend = 2)
  text(0.72, 0.3, "Predicted (+/- Cred. ints.)", adj = 0, cex = 1.25)
  
  # Y-axis labels:
  for (ii in 1:2) {
    par(mar = mar * c(1, 0, 1, 0))
    plot(c(0, 1), c(0, 1), col = NA, xlab = "", ylab = "", axes = FALSE)
    text(0.5, 0.5, c("Fertility", "Number of offspring")[ii], cex = 2, srt = 90)
  }
  
  # X-axis:
  for (ii in 1:2) {
    par(mar = mar * c(0, 1, 0, 1))
    plot(xlim, c(0, 1), col = NA, xlab = "", ylab = "", axes = FALSE)
    Axis(xlim, side = 1, lwd = NA, pos = 1.8, cex.axis = 1.25)
  }
  
  # Main plots:
  idsort <- x$idsort$DIC
  if (sortBy == "PredLoss") {
    idsort <- x$idsort$PredLoss
  }
  for (imod in idsort) {
    for (jj in 1:2) {
      par(mar = mar * c(0, 1, 0, 1))
      plot(c(0, 1), c(0, 1), col = NA, xlab = "", ylab = "", axes = FALSE)
      if (jj == 1) {
        text(0.1, 0.15, modnames[imod], cex = 1.5, adj = 0, xpd = NA)
      }
    }
    par(mar = mar) 
    # Plot fertility:
    dat <- x$runs[[imod]]$data
    alpha <- dat$alpha
    xv <- x$runs[[imod]]$x
    idplf <- which(xv <= max(dat$x))
    fertQuant <- x$runs[[imod]]$fert[idplf, ]
    xv <- xv[idplf]
    plot(xlim, ylim$fert, col = NA, xlab = "", ylab = "", axes = FALSE)
    polygon(dat$alpha + c(xv, rev(xv)), 
            c(fertQuant[, "Lower"], rev(fertQuant[, "Upper"])),
            col = "orange", border = NA)
    lines(xv + alpha, fertQuant[, "Mean"], col = 'red', lwd = 1)
    lines(dat$data$Age[dat$idages] + 0.5, dat$data$Fertility[dat$idages], type = 'b',
          lwd = 2)
    Axis(x = ylim$fert, side = 2, pos = 0, las = 2, cex.axis = 1.25)
    Axis(x = xlim, side = 1, pos = 0, labels = NA)
    for (jj in 1:2) {
      lines(xlim, rep(ylim$fert[jj], 2))
      lines(rep(xlim[jj], 2), ylim$fert)
    }
    
    # Plot predicted number of offspring:
    predQuant <- x$runs[[imod]]$pred
    plot(xlim, ylim$pred, col = NA, xlab = "", ylab = "", axes = FALSE)
    polygon(dat$alpha + c(dat$x, rev(dat$x)), 
            c(predQuant[, "Lower"], rev(predQuant[, "Upper"])),
            col = "orange", border = NA)
    lines(dat$x + alpha, predQuant[, "Mean"], col = 'red', lwd = 1)
    lines(dat$x + alpha, dat$data$nOffspring[dat$idages], col = 1, 
          lwd = 2, type = 'b')
    Axis(x = ylim$pred, side = 2, pos = 0, las = 2, cex.axis = 1.25)
    Axis(x = xlim, side = 1, pos = 0, labels = NA)
    for (jj in 1:2) {
      lines(xlim, rep(ylim$pred[jj], 2))
      lines(rep(xlim[jj], 2), ylim$pred)
    }
    text(xlim[2] - diff(xlim) * 0.35, ylim$pred[2] - diff(ylim$pred) * 0.2,
         sprintf("PL = %s",
                 format(round(x$runs[[imod]]$PredLoss[1,"Deviance"]),
                        big.mark = ",")),
         adj = 0, cex = 1.2)
    text(xlim[2] - diff(xlim) * 0.35, ylim$pred[2] - diff(ylim$pred) * 0.4,
         sprintf("DIC = %s",
                 format(round(x$runs[[imod]]$DIC["DIC"]),
                        big.mark = ",")),
         adj = 0, cex = 1.2)
    
  }
}

# Print for multi-BaFTA:
print.multibafta <- function(x, ...) {
  args <- list(...)
  argNames <- names(args)
  
  if ("digits" %in% argNames) {
    digits <- args$digits
  } else {
    digits <- 3
  }
  cat("Sorted by DIC:\n")
  print(x$gof[x$idsort$DIC, ], digits = digits)
  
  cat("\nSorted by PredLoss:\n")
  print(x$gof[x$idsort$PredLoss, ], digits = digits)
}

# Empty plot:
plot0 <- function(x = c(0, 1), y = c(0, 1), ...) {
  plot(x, y, col = NA, xlab = "", ylab = "", axes = FALSE)
}

# Plot of combined traces and posterior densities:
PlotDensTrace <- function(out, plotDir, figName, saveResults = FALSE,
                          simPars = NULL) {
  # Number of MCMC chains:
  nsim <- out$settings$nsim
  
  # Number of parameters:
  p <- nrow(out$coefficients)
  
  # Parameter names:
  pnames <- rownames(out$coefficients)
  
  # Fancy names:
  labpnames <- c(b0 = expression(theta[0]), b1 = expression(theta[1]), 
                 b2 = expression(theta[2]), b3 = expression(theta[3]),
                 b4 = expression(theta[4]), b5 = expression(theta[5]), 
                 b6 = expression(theta[6]), gamma = expression(gamma),
                 eta = expression(eta), kappa = expression(kappa), 
                 vSd = expression(italic(sigma[v])))
  labpnames <- labpnames[which(names(labpnames) %in% pnames)]
  
  # Converged parameters:
  idSamp <- out$params$idSamp
  if (out$settings$dataType == "indivExtended") {
    idSamp <- c(idSamp, max(idSamp) + 1)
  }
  
  # Converged parameter matrix:
  parsMat <- out$theta[, pnames]
  
  # data for trace plots:
  keep <- seq(1, out$settings$niter, out$settings$thinning)
  nkeep <- length(keep)
  fullPars <- array(0, dim = c(nkeep, p, nsim))
  for (ii in 1:nsim) {
    fullPars[,, ii] <- out$runs[[ii]]$theta[keep, pnames]
  }
  
  # Y limits for traces:
  ylimtr <- sapply(1:p, function(ip) {
    range(fullPars[, ip, ])
  })
  dimnames(ylimtr) <- list(c("low", "upp"), pnames)
  
  # x limits for densities:
  xlimden <- apply(parsMat, 2, function(pv) {
    quantile(pv, c(0.001, 0.999))
  })
  dimnames(xlimden) <- list(c("low", "upp"), pnames)
  
  # Prepare layout:
  ppmat <- matrix(1:(p * 2) + 4, p, 2)
  ppmat <- cbind(ppmat[, 1], rep(2, p), ppmat[, 2])
  laymat <- cbind(c(rep(1, p), 0), rbind(ppmat, c(3, 0, 4)))
  widths <- c(0.2, 1, 0.2, 1)
  heights <- c(rep(0.5, p), 0.2)
  
  # PDF settings:
  whr <- sum(widths) / sum(heights)
  pdfw <- 6
  pdfh <- pdfw / whr
  
  # Margins:
  mar <- c(2, 2, 1, 1)
  
  # Trace colors:
  coltr <- brewer.pal(9, "Set1")[-6][1:nsim]
  
  # Produce plot:
  if (saveResults) {
    pdf(file = sprintf("%s%s.pdf", plotDir, figName), width = pdfw, 
        height = pdfh)
  }
  
  layout(mat = laymat, widths = widths, heights = heights)
  
  # Y label 1:
  par(mar = mar * c(1, 0, 1, 0))
  plot0()
  text(0.5, 0.5, "Parameter posterior density", srt = 90, cex = 2)
  
  # Y label 2:
  par(mar = mar * c(1, 0, 1, 0))
  plot0()
  text(0.5, 0.5, "Parameter traces", srt = 90, cex = 2)
  
  # X label 1:
  par(mar = mar * c(0, 1, 0, 1))
  plot0()
  text(0.5, 0.5, "Parameter value", cex = 2)
  
  # Y label 2:
  par(mar = mar * c(0, 1, 0, 1))
  plot0()
  text(0.5, 0.5, "MCMC step", cex = 2)
  
  # Density plots
  for (ip in 1:p) {
    # Posterior density:
    thev <- parsMat[, ip]
    idth <- which(thev >= xlimden[1, ip] & thev <= xlimden[2, ip])
    thev <- thev[idth]
    thest <- out$coefficients[ip, c(1:4)]
    dthe <- density(thev)
    ylim <- c(0, max(dthe$y))
    idin <- which(dthe$x >= thest["Lower"] & dthe$x <= thest["Upper"])
    idmean <- which(abs(dthe$x - thest["Mean"]) == 
                      min(abs(dthe$x - thest["Mean"])))[1]
    par(mar = mar, las = 1)
    plot(xlimden[, ip], ylim, col = NA, lwd = 2, xlab = "", ylab = "")
    polygon(dthe$x[c(idin, rev(idin))], c(dthe$y[idin], rep(0, length(idin))),
            col = "orange", border = NA)
    lines(dthe$x[rep(idmean, 2)], c(0, dthe$y[idmean]), col = 'white', lwd = 2)
    lines(dthe$x, dthe$y, col = 'orange', lwd = 2)
    if (!is.null(simPars)) {
      lines(x = rep(simPars[ip], 2), y = ylim, col = 'dark red', 
            lwd = 2)
    }
    text(xlimden[1, ip] + diff(xlimden[, ip]) * 0.025, 
         ylim[2] - diff(ylim) * 0.2,
         labpnames[ip], cex = 2, adj = 0)
  }
  
  # Trace plots:
  xlim <- range(keep)
  for (ip in 1:p) {
    par(mar = mar, las = 1)
    plot(xlim, ylimtr[, ip], col = NA, xlab = "", ylab = "")
    for (isim in 1:nsim) {
      lines(keep, fullPars[, ip, isim], col = coltr[isim])
    }
  }
  if (saveResults) dev.off()
  
}
