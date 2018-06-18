import           HNum.Vector                    ( matrix )
import           HNum.CSV
import           Data.List
import           Control.Parallel
import           Control.Concurrent.ParallelIO

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

------------------------------------------
-- Iteration
------------------------------------------
-- Prepare Iteration
type Func = Xi -> Double -> Double

type Xi = Double

data Bosons = Bosons { fh :: Func
                     , fw :: Func
                     , fz :: Func}

data Fermions = Fermions { fu :: Func
                         , fd :: Func
                         , fc :: Func
                         , fs :: Func
                         , ft :: Func
                         , fb :: Func
                         , fe :: Func
                         , fmu :: Func
                         , ftau :: Func
                         , fnu :: Func}

simulExecB :: Xi -> Bosons -> [Double] -> [Double]
simulExecB xi bs xs = gs
 where
  gs = map (apl bs) xs
  apl cs x =
    let hs = fh cs xi x
        ws = fw cs xi x
        zs = fz cs xi x
    in  hs `par` ws `par` zs `pseq` (hs + ws + zs)

simulExecF :: Xi -> Fermions -> [Double] -> [Double]
simulExecF xi ferms xs = gs
 where
  gs = map (apl ferms) xs
  apl hs x =
    let us   = fu hs xi x
        ds   = fd hs xi x
        cs   = fc hs xi x
        ss   = fs hs xi x
        ts   = ft hs xi x
        bs   = fb hs xi x
        es   = fe hs xi x
        mus  = fmu hs xi x
        taus = ftau hs xi x
        nus  = 3 * fnu hs xi x
    in  us
        `par`  ds
        `par`  cs
        `par`  ss
        `par`  ts
        `par`  bs
        `par`  es
        `par`  mus
        `par`  taus
        `par`  nus
        `pseq` (us + ds + cs + ss + ts + bs + es + mus + taus + nus)

zeroF :: Func
zeroF _ _ = 0

type FuncV = Xi -> [Double] -> [Double]

zeroV :: FuncV
zeroV _ xs = replicate (length xs) 0

data Totals = Totals { t1 :: FuncV
                     , t2 :: FuncV
                     , t3 :: FuncV
                     , t4 :: FuncV}

simulExecV :: Xi -> Totals -> [Double] -> [Double]
simulExecV xi = apl
 where
  apl ps ys =
    let x1   = t1 ps xi ys
        x2   = t2 ps xi ys
        x3   = t3 ps xi ys
        x4   = t4 ps xi ys
        x12  = zipWith (+) x1 x2
        x34  = zipWith (+) x3 x4
        xtot = zipWith (+) x12 x34
    in  x1 `par` x2 `par` x3 `par` x4 `pseq` x12 `par` x34 `pseq` xtot

------------------------------------------
-- Two Body
------------------------------------------

-- Decay Width
gphh :: Double -> Double -> Double
gphh xi mPhi
  | mPhi <= 2 * mH
  = 0
  | otherwise
  = xH
    `seq` fac
    `seq` former
    `seq` latter
    `seq` fac
    *     mPhi
    ^     3
    *     former
    *     latter
 where
  xH     = (mH / mPhi) ^ 2
  fac    = factor xi 32
  former = (1 + 2 * xH) ^ 2
  latter = sqrt (1 - 4 * xH)

gpZZ :: Double -> Double -> Double
gpZZ xi mPhi
  | mPhi <= 2 * mZ
  = 0
  | otherwise
  = xZ
    `seq` fac
    `seq` former
    `seq` latter
    `seq` fac
    *     mPhi
    ^     3
    *     former
    *     latter
 where
  xZ     = (mZ / mPhi) ^ 2
  fac    = factor xi 32
  former = 1 - 4 * xZ + 12 * xZ ^ 2
  latter = sqrt (1 - 4 * xZ)

gpWW :: Double -> Double -> Double
gpWW xi mPhi
  | mPhi <= 2 * mW
  = 0
  | otherwise
  = xW
    `seq` fac
    `seq` former
    `seq` latter
    `seq` fac
    *     mPhi
    ^     3
    *     former
    *     latter
 where
  xW     = (mW / mPhi) ^ 2
  fac    = factor xi 16
  former = (1 - 4 * xW + 12 * xW ^ 2) ^ 2
  latter = sqrt (1 - 4 * xW)

-- For Fermion
gpff :: Double -> Double -> Xi -> Double -> Double
gpff mF nC xi mPhi
  | mPhi <= 2 * mF
  = 0
  | otherwise
  = xF `seq` fac `seq` middle `seq` fac * nC * mPhi ^ 3 * xF * middle
 where
  xF     = (mF / mPhi) ^ 2
  fac    = factor xi 8
  middle = sqrt (1 - 4 * xF) ^ 3

-- Fermions
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

-- Bosonic Vectorize
gpBB' :: Xi -> [Double] -> [Double]
gpBB' xi = simulExecB xi bs where bs = Bosons gphh gpWW gpZZ

-- Fermionic Vectorize
gpFF' :: Xi -> [Double] -> [Double]
gpFF' xi = simulExecF xi ferms
  where ferms = Fermions gpUU gpDD gpCC gpSS gpTT gpBB gpEE gpMu gpTau gpNu

-- Two Body Total
gpTot2B' :: Xi -> [Double] -> [Double]
gpTot2B' xi = simulExecV xi ts where ts = Totals gpBB' gpFF' zeroV zeroV

------------------------------------------
-- Three Body
------------------------------------------
-- Three Body
gpqqg :: Double -> Double -> Double
gpqqg xi mPhi = as `seq` fac `seq` 6 * as * fac * mPhi ^ 3
 where
  as  = aS mPhi
  fac = factor xi 4

gpffW :: Double -> Double -> Double -> Double
gpffW nC xi mPhi = fac `seq` gF * nC * fac * mPhi ^ 5
  where fac = 3 * 3 / (4 * sqrt 2) * factor xi (64 * pi)

gpffZ :: Double -> Double -> Double -> Double -> Double -> Double
gpffZ nC vz az xi mPhi =
  fac `seq` middle `seq` 3 / (2 * sqrt 2) * gF * nC * middle * fac * mPhi ^ 5
 where
  fac    = 3 * 3 / (2 * sqrt 2) * factor xi (64 * pi)
  middle = vz ^ 2 + az ^ 2

gpQQW = gpffW 3
gpLLW = gpffW 1

gpFFW' :: Xi -> [Double] -> [Double]
gpFFW' xi = simulExecB xi bs where bs = Bosons gpqqg gpQQW gpLLW

gpFFW'' :: Xi -> [Double] -> [Double]
gpFFW'' xi = simulExecB xi bs where bs = Bosons gpQQW gpLLW zeroF

gpq1Z =
  let vz = gVZ 0.5 (2 / 3)
      az = gAZ 0.5 (2 / 3)
  in  vz `seq` az `seq` gpffZ 3 vz az

gpq2Z =
  let vz = gVZ (-0.5) (-1 / 3)
      az = gAZ (-0.5) (-1 / 3)
  in  vz `seq` az `seq` gpffZ 3 vz az

gpl1Z =
  let vz = gVZ 0.5 0
      az = gAZ 0.5 0
  in  vz `seq` az `seq` gpffZ 1 vz az

gpl2Z =
  let vz = gVZ (-0.5) (-1)
      az = gAZ (-0.5) (-1)
  in  vz `seq` az `seq` gpffZ 1 vz az

gpFFZ' :: Xi -> [Double] -> [Double]
gpFFZ' xi = simulExecF xi fs
 where
  fs = Fermions gpq1Z gpq2Z zeroF zeroF zeroF zeroF gpl1Z gpl2Z zeroF zeroF

gpTot3B' :: Xi -> [Double] -> [Double]
gpTot3B' xi = simulExecV xi ts where ts = Totals gpFFW' gpFFZ' zeroV zeroV

------------------------------------------
-- Four Body
------------------------------------------
-- Four Body
gpWWhh :: Double -> Double -> Double
gpWWhh xi mPhi = fac `seq` fac * mPhi ^ 7
  where fac = xi ^ 2 / (15 * (8 * pi) ^ 5 * vEW ^ 4 * mpl ^ 2)

gpZZhh :: Double -> Double -> Double
gpZZhh xi mPhi = fac `seq` fac * mPhi ^ 7
  where fac = xi ^ 2 / (30 * (8 * pi) ^ 5 * vEW ^ 4 * mpl ^ 2)

gpTot4B' :: Xi -> [Double] -> [Double]
gpTot4B' xi = simulExecB xi bs where bs = Bosons gpWWhh gpZZhh zeroF

-- Total
gpTot xi = simulExecV xi ts where ts = Totals gpTot2B' gpTot3B' gpTot4B' zeroV


-- Branching Ratio
brpqqg :: Xi -> [Double] -> [Double] -> [Double]
brpqqg xi tots mPhis =
  let a = map (gpqqg xi) mPhis in a `pseq` zipWith (/) a tots

brpFF :: Xi -> [Double] -> [Double] -> [Double]
brpFF xi tots mPhis = let a = gpFF' xi mPhis in a `pseq` zipWith (/) a tots

brphh :: Xi -> [Double] -> [Double] -> [Double]
brphh xi tots mPhis =
  let a = map (gphh xi) mPhis in a `pseq` zipWith (/) a tots

brpWWZZ :: Xi -> [Double] -> [Double] -> [Double]
brpWWZZ xi tots mPhis =
  let a = map (apl (gpWW, gpZZ)) mPhis
  in  a `par` tots `pseq` zipWith (/) a tots
 where
  apl (w, z) m =
    let ws = w xi m
        zs = z xi m
    in  ws `par` zs `pseq` ws + zs

brpFFWFFZ :: Xi -> [Double] -> [Double] -> [Double]
brpFFWFFZ xi tots mPhis =
  let a = gpFFW'' xi mPhis in a `pseq` zipWith (/) a tots

brpWWhhZZhh :: Xi -> [Double] -> [Double] -> [Double]
brpWWhhZZhh xi tots mPhis =
  let a = gpTot4B' xi mPhis in a `pseq` zipWith (/) a tots

main
  = let
      v            = [1 .. 1e+6]
      tots         = gpTot 1 [1 .. 1e+6]
      brpqqg'      = map show $ brpqqg 1 tots v
      brpFF'       = map show $ brpFF 1 tots v
      brphh'       = map show $ brphh 1 tots v
      brpWWZZ'     = map show $ brpWWZZ 1 tots v
      brpFFWFFZ'   = map show $ brpFFWFFZ 1 tots v
      brpWWhhZZhh' = map show $ brpWWhhZZhh 1 tots v
    in
      v
      `pseq` tots
      `pseq` brpqqg'
      `par`  brpFF'
      `par`  brphh'
      `par`  brpWWZZ'
      `par`  brpFFWFFZ'
      `par`  brpWWhhZZhh'
      `pseq` do
               parallel_
                 [ writeFile "CSV/total.csv" (intercalate "\n" (map show tots))
                 , writeFile "CSV/brpqqg.csv"    (intercalate "\n" brpqqg')
                 , writeFile "CSV/brpFF.csv"     (intercalate "\n" brpFF')
                 , writeFile "CSV/brphh.csv"     (intercalate "\n" brphh')
                 , writeFile "CSV/brpWWZZ.csv"   (intercalate "\n" brpWWZZ')
                 , writeFile "CSV/brpFFWFFZ.csv" (intercalate "\n" brpFFWFFZ')
                 , writeFile "CSV/brpWWhhZZhh.csv"
                             (intercalate "\n" brpWWhhZZhh')
                 ]
               stopGlobalPool


