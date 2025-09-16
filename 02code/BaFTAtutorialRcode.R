# ============================== CODE METADATA =============================== #
# TITLE: BaFTA tutorial code.
# AUTHOR: Fernando Colchero
# DATE CREATED: 2025-09-10
# DESCRIPTION: Runs BaFTA on aggregated and individual level data with variable
#              IBIs.
# MODEL INCLUDED IN BaFTA: quadratic, PeristeraKostaki, ColcheroMuller, 
#             Hadwiger, gamma, beta, skewNormal, gammaMixture, HadwigerMixture, 
#             skewSymmetric, skewLogistic.
# ================================ CODE START ================================ #
# ================= #
# ==== SET UP: ====
# ================= #
# Installed packages:
instPacks <- installed.packages()[, 1]

# Libraries:
if (!"snowfall" %in% instPacks) {
  install.packages("snowfall")
}

if (!"RColorBrewer" %in% instPacks) {
  install.packages("RColorBrewer")
}

if (!"BaFTA" %in% instPacks) {
  if (!"devtools" %in% instPacks) {
    install.packages("devtools")
  }
  library(devtools)
  
  install_git("https://github.com/fercol/BaFTA", subdir = "pkg/")
}

# Load libraries:
library(BaFTA)
library(RColorBrewer)

# Set working directory (change accordingly):
setwd("~/FERNANDO/MEETINGS/WORKSHOPS/2025.Oxford.BaSTA.BaFTA/BaSTA_BaFTA_Workshop/")

# Load additional functions:
source("02code/BaFTAextraFunctions.R")

# Load aggregated data:
datAggr <- read.csv(file = "03data/BaFTA/speciesFert.csv")

# Load individual data:
datIndivSimp <- read.csv(file = "03data/BaFTA/indivSimpFert.csv")

# Load individual data:
datIndivExt <- read.csv(file = "03data/BaFTA/indivExtFert.csv")

# ====================================== #
# ==== ANALYSIS ON AGGREGATED DATA: ====
# ====================================== #
# Extract species:
species <- unique(datAggr$species)

# Select one of the species:
isp <- 1

# Subset data:
spDat <- datAggr[which(datAggr$species == species[isp]), ]

# ---------------------------- #
# ---- Single model runs: ----
# ---------------------------- #
# Quadratic:
outAggr0 <- bafta(object = spDat, model = "quadratic", nsim = 6, ncpus = 6)
outAggr0
plot(outAggr0)
plot(outAggr0, type = "fertility")
plot(outAggr0, type = "predictive")

# Gamma:
outAggr1 <- bafta(object = spDat, model = "gamma", nsim = 6, ncpus = 6)
outAggr1
plot(outAggr1)
plot(outAggr1, type = "fertility")
plot(outAggr1, type = "predictive")

# ------------------------------ #
# ---- Multiple model runs: ----
# ------------------------------ #
# Models:
models <- c("quadratic", "PeristeraKostaki", "gamma", "beta", "gammaMixture",
            "skewNormal")

# Run multi-BaFTA:
multiout <- multibafta(object = spDat, models = models)

# Print output:
multiout

# Plot output:
plot(multiout, sortBy = "PredLoss")

# ====================================================== #
# ==== INDIVIDUAL LEVEL DATA SEASONAL REPRODUCTION: ====
# ====================================================== #
# -------------------------- #
# ---- Quadratic model: ----
# -------------------------- #
# run analysis:
outIndSimp0 <- bafta(object = datIndivSimp, dataType = "indivSimple", 
                 model = "quadratic", niter = 20000, burnin = 1001, nsim = 4, 
                 ncpus = 4)

# print outIndput to the screen:
outIndSimp0

# plot traces:
plot(outIndSimp0)

# Plot parameter posterior densities:
plot(outIndSimp0, type = "density")

# plot fertility:
plot(outIndSimp0, type = "fertility")

# plot predicted number of offspring:
plot(outIndSimp0, type = "predictive")

# ---------------------- #
# ---- Gamma model: ----
# ---------------------- #
# Run analysis:
outIndSimp1 <- bafta(object = datIndivSimp, dataType = "indivSimple",
                      model = "gamma", niter = 20000, burnin = 1001, nsim = 4, 
                      ncpus = 4)

# print outIndput to the screen:
outIndSimp1

# plot traces:
plot(outIndSimp1)

# Plot parameter posterior densities:
plot(outIndSimp1, type = "density")

# plot fertility:
plot(outIndSimp1, type = "fertility")

# plot predicted number of offspring:
plot(outIndSimp1, type = "predictive")

# Combined plot of posterior densities and traces:
PlotDensTrace(out = outIndSimp1)

# =================================================== #
# ==== INDIVIDUAL LEVEL DATA WITH VARIABLE IBIS: ====
# =================================================== #
# -------------------------- #
# ---- Quadratic model: ----
# -------------------------- #
# run analysis:
outIndExt0 <- bafta(object = datIndivExt, dataType = "indivExtended", 
              model = "quadratic", niter = 20000, burnin = 1001, nsim = 4, 
              ncpus = 4)

# print outIndput to the screen:
outIndExt0

# plot traces:
plot(outIndExt0)

# Plot parameter posterior densities:
plot(outIndExt0, type = "density")

# plot fertility:
plot(outIndExt0, type = "fertility")

# plot predicted number of offspring:
plot(outIndExt0, type = "predictive")

# ---------------------- #
# ---- Gamma model: ----
# ---------------------- #
# Run analysis:
outIndExt1 <- bafta(object = datIndivExt, dataType = "indivExtended", 
                    model = "gamma", niter = 20000, burnin = 1001, nsim = 4, 
                    ncpus = 4)

# print outIndput to the screen:
outIndExt1

# plot traces:
plot(outIndExt1)

# Plot parameter posterior densities:
plot(outIndExt1, type = "density")

# plot fertility:
plot(outIndExt1, type = "fertility")

# plot predicted number of offspring:
plot(outIndExt1, type = "predictive")

# Combined plot of posterior densities and traces:
PlotDensTrace(out = outIndExt1)

# ---------------------------------- #
# ---- Peristera-Kostaki model: ----
# ---------------------------------- #
# Run analysis:
outIndExt2 <- bafta(object = datIndivExt, dataType = "indivExtended", 
              model = "PeristeraKostaki",
              niter = 20000, burnin = 1001, nsim = 4, ncpus = 4)

# print outIndput to the screen:
outIndExt2

# plot traces:
plot(outIndExt2)

# Plot parameter posterior densities:
plot(outIndExt2, type = "density")

# plot fertility:
plot(outIndExt2, type = "fertility")

# plot predicted number of offspring:
plot(outIndExt2, type = "predictive")

# Combined plot of posterior densities and traces:
PlotDensTrace(out = outIndExt2)
