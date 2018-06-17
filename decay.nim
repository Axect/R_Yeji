import math, strutils, sequtils, future

# Mass of Particles
const
  mH = 125.09
  mZ = 91.1876
  mW = 80.385
  mU = 2.2e-3
  mD = 4.7e-3
  mS = 9.6e-2
  mC = 1.28
  mT = 173.1
  mB = 4.18
  mE = 5.109989461e-4
  mMu = 1.056583745e-5
  mTau = 1.77686
  mNu = 0
  M = 2.435e18
  vEW = 246
  GF = 1.1663787e-5
  sw = sqrt(0.23)
  cw = sqrt(1-sw^2)
  aW = 1/28
  g2 = sqrt(4*PI*aW)
  GeVinSecond = 1.52e+24

proc aSInv(muS: float64): float64 =
  result = 1/0.118 - (2*4*PI)/(16*PI^2) 
  
