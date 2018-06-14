import           HNum.Vector
import           HNum.Stats
import           HNum.CSV
import           Data.List

-- Consts
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
mpl = 2.435e18
vEW = 246
gF = 1.1663787e-5
sw = sqrt 0.23
cw = sqrt (1 - sw ^ 2)
aW = 1 / 28
g2 = sqrt (4 * pi * aW)
geVinSec = 1.52e+24

aSInv :: Double -> Double
aSInv muS =
  1 / 0.118 - (2 * 4 * pi) / (16 * pi ^ 2) * (-11 + 4 / 3 * 6) * log (muS / 91)

aS :: Double -> Double
aS mPhi = 1 / aSInv mPhi

gVZ :: Double -> Double -> Double
gVZ tL3 q = g2 / (2 * cw) * (tL3 - 2 * q * sw ^ 2)

gAZ :: Double -> Double -> Double
gAZ tL3 _ = g2 / (2 * cw) * tL3

factor :: Double -> Double -> Double
factor xi num = xi ^ 2 / (num * pi * mpl ^ 2)

simulApply2 :: [Double -> Double] -> Double -> Double
simulApply2 xs d = foldr (\x -> (+) (x d)) 0 xs

simulApply3 :: [Double -> Double -> Double] -> Double -> Double -> Double
simulApply3 xs d1 d2 = foldr (\x -> (+) (x d1 d2)) 0 xs
------------------------------------------
-- Two Body
------------------------------------------

-- Decay Width
gphh :: Double -> Double -> Double
gphh xi mPhi
  | mPhi <= 2 * mH = 0
  | otherwise = factor xi 32 * mPhi ^ 3 * (1 + 2 * xH) ^ 2 * sqrt (1 - 4 * xH)
  where xH = (mH / mPhi) ^ 2

gpZZ :: Double -> Double -> Double
gpZZ xi mPhi
  | mPhi <= 2 * mZ = 0
  | otherwise = factor xi 32 * mPhi ^ 3 * (1 - 4 * xZ + 12 * xZ ^ 2) * sqrt
    (1 - 4 * xZ)
  where xZ = (mZ / mPhi) ^ 2

gpWW :: Double -> Double -> Double
gpWW xi mPhi
  | mPhi <= 2 * mW = 0
  | otherwise = factor xi 16 * mPhi ^ 3 * (1 - 4 * xW + 12 * xW ^ 2) ^ 2 * sqrt
    (1 - 4 * xW)
  where xW = (mW / mPhi) ^ 2

-- For Fermion
gpff :: Double -> Double -> Double -> Double -> Double
gpff mF nC xi mPhi
  | mPhi <= 2 * mF = 0
  | otherwise      = factor xi 8 * nC * mPhi ^ 3 * xF * sqrt (1 - 4 * xF) ^ 3
  where xF = (mF / mPhi) ^ 2

gpUU = gpff mU 3
gpDD = gpff mD 3
gpCC = gpff mC 3
gpSS = gpff mS 3
gpTT = gpff mT 3
gpBB = gpff mB 3
gpEE = gpff mE 1
gpMu = gpff mMu 1
gpTau = gpff mTau 1
gpNu = gpff mNu 1

gpFF :: Double -> Double -> Double
gpFF =
  simulApply3 [gpUU, gpDD, gpCC, gpSS, gpTT, gpBB, gpEE, gpMu, gpTau, gpNu]

gpTot2B :: Double -> Double -> Double
gpTot2B = simulApply3 [gphh, gpZZ, gpWW, gpFF]

-- Three Body
gpqqg :: Double -> Double -> Double
gpqqg xi mPhi = 6 * aS mPhi * factor xi 4 * mPhi ^ 3

gpffW :: Double -> Double -> Double -> Double
gpffW nC xi mPhi =
  3 * 3 / (4 * sqrt 2) * gF * nC * factor xi (64 * pi) * mPhi ^ 5

gpffZ :: Double -> Double -> Double -> Double -> Double -> Double
gpffZ nC vz az xi mPhi =
  3
    * 3
    / (2 * sqrt 2)
    * gF
    * nC
    * (vz ^ 2 + az ^ 2)
    * factor xi (64 * pi)
    * mPhi
    ^ 5

gpqqW = gpffW 3
gpllW = gpffW 1
gpFFW = simulApply3 [gpqqW, gpllW]
gpFFZ = simulApply3
  [ gpffZ 3 (gVZ 0.5 (2 / 3))       (gAZ 0.5 (2 / 3))
  , gpffZ 3 (gVZ (-1 / 2) (-1 / 3)) (gAZ (-1 / 2) (-1 / 3))
  , gpffZ 1 (gVZ 0.5 0)             (gAZ 0.5 0)
  , gpffZ 1 (gVZ (-1 / 2) (-1))     (gAZ (-1 / 2) (-1))
  ]

gpTot3B = simulApply3 [gpqqg, gpFFW, gpFFZ]

-- Four Body
gpWWhh :: Double -> Double -> Double
gpWWhh xi mPhi = xi ^ 2 / (15 * (8 * pi) ^ 6 * vEW ^ 4 * mpl ^ 2) * mPhi ^ 7

gpZZhh :: Double -> Double -> Double
gpZZhh xi mPhi = xi ^ 2 / (30 * (8 * pi) ^ 5 * vEW ^ 4 * mpl ^ 2) * mPhi ^ 7

gpTot4B = simulApply3 [gpWWhh, gpZZhh]

-- Total
gpTot = simulApply3 [gpTot2B, gpTot3B, gpTot4B]


main = do
  let
    !v           = vec [1 .. 1e6]
    !vTot        = gpTot 1 <$> v
    !vpqqg       = gpqqg 1 <$> v
    !vpFF        = gpFF 1 <$> v
    !vpWW        = gpWW 1 <$> v
    !vpZZ        = gpZZ 1 <$> v
    !vphh        = gphh 1 <$> v
    !vpFFW       = gpFFW 1 <$> v
    !vpFFZ       = gpFFZ 1 <$> v
    !vpWWhh      = gpWWhh 1 <$> v
    !vpZZhh      = gpZZhh 1 <$> v
    !brpqqg      = toList $ vpqqg / vTot
    !brpFF       = toList $ vpFF / vTot
    !brpWWZZ     = toList $ (vpZZ + vpWW) / vTot
    !brphh       = toList $ vphh / vTot
    !brpFFWFFZ   = toList $ (vpFFW + vpFFZ) / vTot
    !brpWWhhZZhh = toList $ (vpWWhh + vpZZhh) / vTot
    !df          = DataFrame
      { header = [ "Total"
                 , "brpqqg"
                 , "brpFF"
                 , "brpWWZZ"
                 , "brphh"
                 , "brpFFWFFZ"
                 , "brpWWhhZZhh"
                 ]
      , dat    = matrix
        [toList vpqqg, brpqqg, brpFF, brpWWZZ, brphh, brpFFWFFZ, brpWWhhZZhh]
      }
  writeCSV "test.csv" df

