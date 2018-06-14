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
aS <- function (mEta) 1 / (aSInv(mEta))
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
Gehh <- function (mEta, xi) {
  if (mEta <= 2*mH) {
    return(0)
  } else {
    xH <- (mH / mEta)^2
    return(xi^2/(32*pi*vEW^2) * mEta^3 * (1+2*xH)^2 * sqrt(1 - 4*xH))
  }
}

GeZZ <- function (mEta, xi) {
  if (mEta <= 2*mZ) {
    return(0)
  } else {
    xZ <- (mZ / mEta)^2
    return(xi^2/(32*pi*vEW^2) * mEta^3 * (1 - 4*xZ + 12*xZ^2) * sqrt(1 - 4*xZ))
  }
}

GeWW <- function (mEta, xi) {
  if (mEta <= 2*mW) {
    return(0)
  } else {
    xW <- (mW / mEta)^2
    return(xi^2/(16*pi*vEW^2) * mEta^3 * (1 - 4*xW + 12*xW^2) * sqrt(1 - 4*xW))
  }
}

Geff <- function (mEta, xi, mF, NC) {
  if (mEta <= 2*mF) {
    return(0)
  } else {
    xF <- (mF / mEta)^2
    return(NC * xi^2/(8*pi*vEW^2) * mEta^3 * xF * sqrt(1 - 4*xF)^3)
  }
}

GeU <- function (mEta, xi) Geff(mEta, xi, mU, 3)
GeD <- function (mEta, xi) Geff(mEta, xi, mD, 3)
GeC <- function (mEta, xi) Geff(mEta, xi, mC, 3)
GeS <- function (mEta, xi) Geff(mEta, xi, mS, 3)
GeT <- function (mEta, xi) Geff(mEta, xi, mT, 3)
GeB <- function (mEta, xi) Geff(mEta, xi, mB, 3)
GeE <- function (mEta, xi) Geff(mEta, xi, mE, 1)
GeMu <- function (mEta, xi) Geff(mEta, xi, mMu, 1)
GeTau <- function (mEta, xi) Geff(mEta, xi, mTau, 1)
GeNu <- function (mEta, xi) 3*Geff(mEta, xi, mNu, 1)

# Vectorize Functions
gehh <- Vectorize(Gehh)
geZZ <- Vectorize(GeZZ)
geWW <- Vectorize(GeWW)
geU <- Vectorize(GeU)
geD <- Vectorize(GeD)
geC <- Vectorize(GeC)
geS <- Vectorize(GeS)
geT <- Vectorize(GeT)
geB <- Vectorize(GeB)
geE <- Vectorize(GeE)
geMu <- Vectorize(GeMu)
geTau <- Vectorize(GeTau)
geNu <- Vectorize(GeNu)
geFF <- function (mEta, xi) {
  return(
    geU(mEta,xi) + geD(mEta,xi)
    + geC(mEta,xi) + geS(mEta,xi)
    + geT(mEta,xi) + geB(mEta,xi)
    + geE(mEta,xi) + geMu(mEta,xi)
    + geTau(mEta,xi) + geNu(mEta,xi)
  )
}
geTot2B <- function (mEta, xi) {
  return(
    gehh(mEta,xi) + geZZ(mEta,xi)
    + geWW(mEta,xi) + geFF(mEta,xi)
  )
}

#------------------------------------------
# Three Body
#------------------------------------------
Geqqg <- function (mEta, xi) {
  return(6 * aS(mEta) * xi^2 / (4*pi^2*vEW^2) * mEta^3)
}

GeffW <- function (mEta, xi, NC) {
  return(3/(4*sqrt(2)) * GF * NC * xi^2 / ((4*pi)^3 * vEW^2) * mEta^5)
}

GeffZ <- function (mEta, xi, NC, vz, az) {
  return(3/(2*sqrt(2)) * GF * NC * (vz^2 + az^2) * xi^2 / ((4*pi)^3 * vEW^2) * mEta^5)
}

GeQQW <- function (mEta, xi) GeffW(mEta, xi, 3)
GellW <- function (mEta, xi) GeffW(mEta, xi, 1)
GeFFZ <- function (mEta, xi) {
  return(
         3*GeffZ(mEta, xi, 3, gVZ(1/2,2/3), gAZ(1/2,2/3))
         + 3*GeffZ(mEta, xi, 3, gVZ(-1/2, -1/3), gAZ(-1/2, -1/3))
         + 3*GeffZ(mEta, xi, 1, gVZ(1/2, 0), gAZ(1/2, 0))
         + 3*GeffZ(mEta, xi, 1, gVZ(-1/2, -1), gAZ(-1/2, -1))
  )
}

# Vectorized
geqqg <- Vectorize(Geqqg)
geQQW <- Vectorize(GeQQW)
gellW <- Vectorize(GellW)
geFFW <- function (mEta, xi) 3*geQQW(mEta, xi) + 3*gellW(mEta, xi)
geFFZ <- Vectorize(GeFFZ)
geTot3B <- function (mEta, xi) {
  return(
    geqqg(mEta, xi) + geFFW(mEta, xi) + geFFZ(mEta, xi)
  )
}


#------------------------------------------
# Four Body
#------------------------------------------
GeffWh <- function (mEta, xi, NC) {
  return(3*sqrt(2)/160 * GF * NC * xi^2/((4*pi)^5 * vEW^4) * mEta^7)
}

GeffZh <- function (mEta, xi, NC, vz, az) {
  return(3*sqrt(2)/80 * GF * NC * (vz^2 + az^2) * xi^2 /((4*pi)^5 * vEW^4) * mEta^7)
}

GeQQWh <- function (mEta, xi) GeffWh(mEta, xi, 3)
GellWh <- function (mEta, xi) GeffWh(mEta, xi, 1)
GeFFZh <- function (mEta, xi) {
  return(
         3*GeffZh(mEta, xi, 3, gVZ(1/2,2/3), gAZ(1/2,2/3))
         + 3*GeffZh(mEta, xi, 3, gVZ(-1/2, -1/3), gAZ(-1/2, -1/3))
         + 3*GeffZh(mEta, xi, 1, gVZ(1/2, 0), gAZ(1/2, 0))
         + 3*GeffZh(mEta, xi, 1, gVZ(-1/2, -1), gAZ(-1/2, -1))
         )
}
# Vectorized
geQQWh <- Vectorize(GeQQWh)
gellWh <- Vectorize(GellWh)
geFFZh <- Vectorize(GeFFZh)
geTot4B <- function (mEta, xi) {
  return(geQQWh(mEta, xi) + gellWh(mEta, xi) + geFFZh(mEta, xi))
}

#------------------------------------------
# Five Body
#------------------------------------------
GeWWhhh <- function (mEta, xi) {
  return(2/(75*(8*pi)^7 * vEW^8) * xi^2 * mEta^9)
}

GeZZhhh <- function (mEta, xi) {
  return(1/(75*(8*pi)^7 * vEW^8) * xi^2 * mEta^9)
}

# Vectorized
geWWhhh <- Vectorize(GeWWhhh)
geZZhhh <- Vectorize(GeZZhhh)
geTot5B <- function (mEta, xi) {
  return(
         geWWhhh(mEta, xi) + geZZhhh(mEta, xi)
         )
}

#------------------------------------------
# Total Width
#------------------------------------------
geTot <- function (mEta, xi) {
  geTot2B(mEta,xi) + geTot3B(mEta, xi) + geTot4B(mEta,xi) + geTot5B(mEta, xi)
}



#------------------------------------------
# Branching ratios
#------------------------------------------
Breqqg <- function (mEta, xi) {
  geqqg(mEta, xi) / geTot(mEta, xi)
}

BreFF <- function (mEta, xi) {
  geFF(mEta, xi) / geTot(mEta, xi)
}

BreWWZZ <- function (mEta, xi) {
  (geZZ(mEta, xi) + geWW(mEta, xi)) / geTot(mEta, xi)
}

Brehh <- function (mEta, xi) {
  gehh(mEta, xi) / geTot(mEta, xi)
}

BreFFWFFZ <- function (mEta, xi) {
  (geFFW(mEta,xi) + geFFZ(mEta, xi)) / geTot(mEta, xi)
}

BreFFWhFFZh <- function (mEta, xi) {
  (geQQWh(mEta, xi) +gellWh(mEta, xi) + geFFZh(mEta, xi)) / geTot(mEta, xi)
}

#BreWWhhZZhh <- function (mEta, xi) {
#  (geWWhh(mEta, xi) + geZZhh(mEta, xi)) / geTot(mEta, xi)
#}

BreWWhhhZZhhh <- function (mEta, xi) {
  (geWWhhh(mEta, xi) + geZZhhh(mEta, xi)) / geTot(mEta, xi)
}
# Insert Value
breqqg <- Breqqg(1:1e+6, 1)
breFF <- BreFF(1:1e+6, 1)
breWWZZ <- BreWWZZ(1:1e+6, 1)
brehh <- Brehh(1:1e+6, 1)
breFFWFFZ <- BreFFWFFZ(1:1e+6, 1)
#breWWhhZZhh <- BreWWhhZZhh(1:1e+6, 1)
breFFWhFFZh <- BreFFWhFFZh(1:1e+6, 1)
breWWhhhZZhhh <- BreWWhhhZZhhh(1:1e+6, 1)

tw1 <- 1 / (geTot(1:1e+6, 1) * GeVinSecond)
tw8 <- 1 / (geTot(1:1e+6, 1e-8) * GeVinSecond)
tw16 <- 1 / (geTot(1:1e+6, 1e-16) * GeVinSecond)

# Plot
png("/home/kavis/Documents/Project/R_Project/Yeji/Figure2.png")
plot(breqqg, log='x', type='l', col='red')
lines(breFF, col='pink')
lines(breWWZZ, col='orange')
lines(brehh, col='blue')
lines(breFFWFFZ, col='green')
lines(breFFWhFFZh, col='purple')
#lines(breWWhhZZhh, col='green')
lines(breWWhhhZZhhh, col='brown')
dev.off()

png("/home/kavis/Documents/Project/R_Project/Yeji/Figure2-2.png")
plot(tw1, log='xy', type='l', col='red')
lines(tw8, col='blue')
lines(tw16, col='green')
dev.off()
