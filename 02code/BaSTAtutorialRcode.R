# ============================== CODE METADATA =============================== #
# TITLE: BaSTA tutorial code.
# AUTHOR: Fernando Colchero
# DATE CREATED: 2025-09-10
# DESCRIPTION: Runs BaSTA on CMR and census data .
# MODEL INCLUDED IN BaSTA: 
#            - Models: Constant mortality (EX), Gompertz (GO), 
#                          Weibull (WE), logistic (LO).
#            - Shapes: simple, Makeham, bathtub
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
  
  install_git("https://github.com/fercol/BaSTA2.0", subdir = "pkg/")
}

# Load libraries:
library(BaSTA)
library(snowfall)
library(RColorBrewer)

# Set working directory (change accordingly):
setwd("~/FERNANDO/MEETINGS/WORKSHOPS/2025.Oxford.BaSTA.BaFTA/BaSTA_BaFTA_Workshop/")

# Load extra function:
source("02code/BaSTAextraFunctions.R")

# ==================== #
# ==== LOAD DATA: ====
# ==================== #
# Sparrowhawk: 
sphDat <- read.csv(file = "03data/BaSTA/SparrowhawkData.csv")

# Elephant data (from Crawley et al (2020) Sci Reps):
eleDat <- read.csv(file = "03data/BaSTA/elephantData.csv")

# =============================== #
# ==== SPARROWHAWK CMR DATA: ====
# =============================== #
# ---------------------------- #
# ---- Single model runs: ----
# ---------------------------- #
# Weibull model:
outSph1 <- basta(object = sphDat, studyStart = 1971, studyEnd = 1999, 
                model = "WE", shape = "simple",
                formulaMort = ~ sex - 1, parallel = TRUE, nsim = 4, ncpus = 4)

outSph1
plot(outSph1)
plot(outSph1, type = 'gof')

# Logistic model:
outSph2 <- basta(object = sphDat, studyStart = 1971, studyEnd = 1999, 
                model = "LO", shape = "simple",
                formulaMort = ~ sex - 1, parallel = TRUE, nsim = 4, ncpus = 4)

outSph2
plot(outSph2)
plot(outSph2, type = 'gof')


# ------------------------------ #
# ---- Multiple model runs: ----
# ------------------------------ #
# Run multiple models:
multiout <- multibasta(object = sphDat, studyStart = 1971, studyEnd = 1999,
                       models = c("GO", "WE", "LO"), shapes = "simple", 
                       formulaMort = ~ sex - 1, 
                       parallel = TRUE, nsim = 4, ncpus = 4)

# Summary of results:
summary(multiout)

# =============================== #
# ==== ELEPHANT CENSUS DATA: ====
# =============================== #
# Check data consistency:
eleCheck <- DataCheck(object = eleDat, dataType = 'census')

# Print the the console results of check:
summary(eleCheck)

# Run test:
outEle1 <- basta(object = eleDat, dataType = "census", model = "GO", 
              shape = "bathtub", formulaMort = ~ Sex - 1, 
              parallel = TRUE, ncpus = 4, nsim = 4)

outEle1
plot(outEle1)
plot(outEle1, type = 'gof')

# Run with minAge = 10:
outEle2 <- basta(object = eleDat, dataType = "census", model = "GO", 
              shape = "bathtub", formulaMort = ~ Sex - 1, minAge = 10,
              parallel = TRUE, ncpus = 4, nsim = 4)

outEle2
plot(outEle2)
plot(outEle2, type = 'gof')
plot(outEle2, type = "demorates")