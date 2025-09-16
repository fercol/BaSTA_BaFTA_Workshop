# ============================== CODE METADATA =============================== #
# AUTHOR: Fernando Colchero
# DATE CREATED: 2025-10-11
# DESCRIPTION: Additional functions for BaSTA analyses.
# NOTES: 
# ================================ CODE START ================================ #
# Create min max columns: 
CalcMinMax <- function(object = object, studyStart, studyEnd) {
  # Study year vector:
  study <- studyStart:studyEnd
  
  # Number of years:
  nStudy <- studyEnd - studyStart + 1
  
  # Maximum age:
  maxAge <- 22
  
  # Birth and death years:
  if (inherits(object, "data.frame")) {
    birth <- object[[2]]
    death <- object[[3]]
  } else {
    birth <- object[, 2]
    death <- object[, 3]
  }
  
  # capture matrix:
  capMat <- object[, 1:nStudy + 3]
  
  # years capture matrix:
  yCapMat <- t(t(capMat) * study)
  
  # First and last capture year:
  flCap <- t(apply(yCapMat, 1, function(icv) {
    maxicv <- study[which(icv == max(icv))]
    minicv <- study[which(icv == min(icv[which(icv > 0)]))]
    return(c(first = minicv, last = maxicv))
  }))
  
  # index of unknown births:
  idub <- which(is.na(birth) & !is.na(death))
  
  # index of unknown deaths:
  idud <- which(is.na(death) & !is.na(birth))
  
  # index of unknown birth death:
  idubd <- which(is.na(birth) & is.na(death))
  
  # Min-max birth:
  minBirth <- maxBirth <- birth
  if (length(idub) > 0) {
    minBirth[idub] <- death[idub] - maxAge
    maxBirth[idub] <- flCap[idub, 1] - 1
  }
  
  # Min-max death:
  minDeath <- maxDeath <- death
  if (length(idud) > 0) {
    maxDeath[idud] <- birth[idud] + maxAge
    minDeath[idud] <- flCap[idud, 2]
  }
  
  # both missing:
  if (length(idubd) > 0) {
    maxBirth[idubd] <-  flCap[idubd, 1] - 1
    minDeath[idubd] <- flCap[idubd, 2]
    minBirth[idubd] <- minDeath[idubd] - maxAge
    maxDeath[idubd] <- maxBirth[idubd] + maxAge
  }
  
  # Correct min-max discrepancies:
  idminmax <- which(minBirth > maxBirth)
  minBirth[idminmax] <- maxBirth[idminmax]
  idminmax <- which(minDeath > maxDeath)
  maxDeath[idminmax] <- minDeath[idminmax]
  
  # Include new columns to object:
  object <- cbind(object, Min.Birth = minBirth, Max.Birth = maxBirth,
                  Min.Death = minDeath, Max.Death = maxDeath)
  return(object)
}
