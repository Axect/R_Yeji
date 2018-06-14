# Mass of Particles
mH <- 125.09
mZ <- 91.1876
mW <- 80.385
mU <- 2.2e-3
mD <- 4.7e-3
mS <- 9.6e-2
mC <- 1.28
mT <- 173.1
mB <- 4.18
mE <- 5.109989461e-4
mMu <- 1.056583745e-5
mTau <- 1.77686
mNu <- 0
M <- 2.435e18
vEW <- 246
GF <- 1.1663787e-5
sw <- sqrt(0.23)
cw <- sqrt(1-sw^2)
aW <- 1/28
g2 <- sqrt(4*pi*aW)
GeVinSecond <- 1.52e+24

# Parameter
aSInv <- function (muS) {
  return(1/0.118 - (2*4*pi)/(16*pi^2) * (-11 + 4/3 * 6) * log(muS / 91))
}
aS <- function (mPhi) 1 / (aSInv(mPhi))
gVZ <- function (TL3, Q) {
  return(g2/(2*cw) * (TL3 - 2*Q*sw^2))
}
gAZ <- function (TL3, Q) {
  return(g2/(2*cw) * TL3)
}

#------------------------------------------
# Two Body
#------------------------------------------
# Decay Width
Gphh <- function (mPhi, xi) {
  if (mPhi <= 2*mH) {
    return(0)
  } else {
    xH <- (mH / mPhi)^2
    return(xi^2/(32*pi*M^2) * mPhi^3 * (1+2*xH)^2 * sqrt(1 - 4*xH))
  }
}

GpZZ <- function (mPhi, xi) {
  if (mPhi <= 2*mZ) {
    return(0)
  } else {
    xZ <- (mZ / mPhi)^2
    return(xi^2/(32*pi*M^2) * mPhi^3 * (1 - 4*xZ + 12*xZ^2) * sqrt(1 - 4*xZ))
  }
}

GpWW <- function (mPhi, xi) {
  if (mPhi <= 2*mW) {
    return(0)
  } else {
    xW <- (mW / mPhi)^2
    return(xi^2/(16*pi*M^2) * mPhi^3 * (1 - 4*xW + 12*xW^2) * sqrt(1 - 4*xW))
  }
}

Gpff <- function (mPhi, xi, mF, NC) {
  if (mPhi <= 2*mF) {
    return(0)
  } else {
    xF <- (mF / mPhi)^2
    return(NC * xi^2/(8*pi*M^2) * mPhi^3 * xF * sqrt(1 - 4*xF)^3)
  }
}

GpU <- function (mPhi, xi) Gpff(mPhi, xi, mU, 3)
GpD <- function (mPhi, xi) Gpff(mPhi, xi, mD, 3)
GpC <- function (mPhi, xi) Gpff(mPhi, xi, mC, 3)
GpS <- function (mPhi, xi) Gpff(mPhi, xi, mS, 3)
GpT <- function (mPhi, xi) Gpff(mPhi, xi, mT, 3)
GpB <- function (mPhi, xi) Gpff(mPhi, xi, mB, 3)
GpE <- function (mPhi, xi) Gpff(mPhi, xi, mE, 1)
GpMu <- function (mPhi, xi) Gpff(mPhi, xi, mMu, 1)
GpTau <- function (mPhi, xi) Gpff(mPhi, xi, mTau, 1)
GpNu <- function (mPhi, xi) 3*Gpff(mPhi, xi, mNu, 1)

# Vectorize Functions
gphh <- Vectorize(Gphh)
gpZZ <- Vectorize(GpZZ)
gpWW <- Vectorize(GpWW)
gpU <- Vectorize(GpU)
gpD <- Vectorize(GpD)
gpC <- Vectorize(GpC)
gpS <- Vectorize(GpS)
gpT <- Vectorize(GpT)
gpB <- Vectorize(GpB)
gpE <- Vectorize(GpE)
gpMu <- Vectorize(GpMu)
gpTau <- Vectorize(GpTau)
gpNu <- Vectorize(GpNu)
gpFF <- function (mPhi, xi) {
  return(
    gpU(mPhi,xi) + gpD(mPhi,xi) 
    + gpC(mPhi,xi) + gpS(mPhi,xi) 
    + gpT(mPhi,xi) + gpB(mPhi,xi)
    + gpE(mPhi,xi) + gpMu(mPhi,xi) 
    + gpTau(mPhi,xi) + gpNu(mPhi,xi) 
  )
} 
gpTot2B <- function (mPhi, xi) {
  return(
    gphh(mPhi,xi) + gpZZ(mPhi,xi)
    + gpWW(mPhi,xi) + gpFF(mPhi,xi)
  )
}

#------------------------------------------
# Three Body
#------------------------------------------
Gpqqg <- function (mPhi, xi) {
  return(6 * aS(mPhi) * xi^2 / (4*pi^2*M^2) * mPhi^3)
}

GpffW <- function (mPhi, xi, NC) {
  return(3/(4*sqrt(2)) * GF * NC * xi^2 / ((4*pi)^3 * M^2) * mPhi^5)
}

GpffZ <- function (mPhi, xi, NC, vz, az) {
  return(3/(2*sqrt(2)) * GF * NC * (vz^2 + az^2) * xi^2 / ((4*pi)^3 * M^2) * mPhi^5)
}

GpQQW <- function (mPhi, xi) GpffW(mPhi, xi, 3)
GpllW <- function (mPhi, xi) GpffW(mPhi, xi, 1)
GpFFZ <- function (mPhi, xi) {
  return(
         3*GpffZ(mPhi, xi, 3, gVZ(1/2,2/3), gAZ(1/2,2/3))
         + 3*GpffZ(mPhi, xi, 3, gVZ(-1/2, -1/3), gAZ(-1/2, -1/3))
         + 3*GpffZ(mPhi, xi, 1, gVZ(1/2, 0), gAZ(1/2, 0))
         + 3*GpffZ(mPhi, xi, 1, gVZ(-1/2, -1), gAZ(-1/2, -1))
  )
}

# Vectorized
gpqqg <- Vectorize(Gpqqg)
gpQQW <- Vectorize(GpQQW)
gpllW <- Vectorize(GpllW)
gpFFW <- function (mPhi, xi) 3*gpQQW(mPhi, xi) + 3*gpllW(mPhi, xi)
gpFFZ <- Vectorize(GpFFZ)
gpTot3B <- function (mPhi, xi) {
  return(
    gpqqg(mPhi, xi) + gpFFW(mPhi, xi) + gpFFZ(mPhi, xi)
  )
}

#------------------------------------------
# Four Body
#------------------------------------------
GpWWhh <- function (mPhi, xi) {
  return(xi^2 / (15*(8*pi)^5*vEW^4*M^2) * mPhi^7)
}

GpZZhh <- function (mPhi, xi) {
  return(xi^2/(30*(8*pi)^5*vEW^4*M^2) * mPhi^7)
}

# Vectorized
gpWWhh <- Vectorize(GpWWhh)
gpZZhh <- Vectorize(GpZZhh)
gpTot4B <- function (mPhi, xi) {
  return(gpWWhh(mPhi, xi) + gpZZhh(mPhi, xi))
}


#------------------------------------------
# Total Width
#------------------------------------------
gpTot <- function (mPhi, xi) {
  gpTot2B(mPhi,xi) + gpTot3B(mPhi, xi) + gpTot4B(mPhi,xi)
}


# Branching Ratio
Brpqqg <- function (mPhi, xi) {
  gpqqg(mPhi,xi) / gpTot(mPhi,xi)
}

BrpFF <- function (mPhi, xi) {
  gpFF(mPhi,xi) / gpTot(mPhi,xi)
}

BrpWWZZ <- function (mPhi, xi) {
  (gpZZ(mPhi, xi) + gpWW(mPhi,xi)) / gpTot(mPhi,xi)
}

Brphh <- function (mPhi, xi) {
  gphh(mPhi,xi) / gpTot(mPhi,xi)
}

BrpFFWFFZ <- function (mPhi, xi) {
  (gpFFW(mPhi,xi) + gpFFZ(mPhi,xi)) / gpTot(mPhi,xi)
}

BrpWWhhZZhh <- function (mPhi, xi) {
  (gpWWhh(mPhi,xi) + gpZZhh(mPhi,xi)) / gpTot(mPhi,xi)
}


# Insert Value
brpqqg <- Brpqqg(1:1e+6, 1)
brpFF <- BrpFF(1:1e+6, 1)
brpWWZZ <- BrpWWZZ(1:1e+6, 1)
brphh <- Brphh(1:1e+6, 1)
brpFFWFFZ <- BrpFFWFFZ(1:1e+6, 1)
brpWWhhZZhh <- BrpWWhhZZhh(1:1e+6, 1)

tw1 <- 1 / (gpTot(1:1e+6, 1) * GeVinSecond)
tw8 <- 1 / (gpTot(1:1e+6, 1e-8) * GeVinSecond)
tw16 <- 1 / (gpTot(1:1e+6, 1e-16) * GeVinSecond)

# Plot
png("/home/kavis/Documents/Project/R_Project/Yeji/Figure1.png")
plot(brpqqg, log='x', type='l', col='red')
lines(brpFF, col='pink')
lines(brpWWZZ, col='orange')
lines(brphh, col='blue')
lines(brpFFWFFZ, col='green')
lines(brpWWhhZZhh, col='brown')
dev.off()

png("Figure1-2.png")
plot(tw1, log='xy', ylim=c(1e+5,1e+45), type='l', col='red')
lines(tw8, col='blue')
lines(tw16, col='green')
dev.off()
