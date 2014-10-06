{-# LANGUAGE MultiParamTypeClasses, FlexibleInstances, UndecidableInstances #-}
{-# LANGUAGE DataKinds, KindSignatures, GADTs #-}
{-# LANGUAGE ScopedTypeVariables #-}

{- | A straightforward LinAlg backend using the hmatrix library. -}

module Numeric.LinAlg.HMatrix (
    module Numeric.LinAlg
  , HArr, hMat, hVec
) where

import Data.Proxy (Proxy (Proxy))

import GHC.TypeLits (SomeNat (SomeNat), someNatVal)

import Numeric.LinAlg
import Numeric.LinAlg.SNat
import qualified Numeric.LinAlg.Vect as V
import Numeric.LinAlg.Vect (Vect (MkVect))
import qualified Numeric.LinearAlgebra as H
import Numeric.LinearAlgebra.HMatrix (Numeric)

import Foreign.Storable (Storable)

data HArr :: Dim -> * -> * where
  HMat :: SNat m -> SNat n -> H.Matrix e -> HArr (M m n) e
  HVec :: SNat n ->           H.Vector e -> HArr (V n  ) e

data SomeMatrix :: * -> * where
  SomeMatrix :: HArr (M m n) e -> SomeMatrix e

data SomeVector :: * -> * where
  SomeVector :: HArr (V n  ) e -> SomeVector e

hMat :: H.Matrix e -> SomeMatrix e
hMat x = let (m, n) = (H.rows x, H.cols x) in 
  case (someNatVal (fromIntegral m), someNatVal (fromIntegral n)) of
    (Just (SomeNat i), Just (SomeNat j)) -> SomeMatrix (HMat (lit i) (lit j) x)

hVec :: Storable e => H.Vector e -> SomeVector e
hVec x = let n = (H.dim x) in 
  case someNatVal (fromIntegral n) of
    Just (SomeNat i) -> SomeVector (HVec (lit i) x)

instance (H.Container H.Vector k) => Eq (HArr (M m n) k) where
  HMat _ _ x == HMat _ _ y = x == y


instance (Eq k, Storable k) => Eq (HArr (V n) k) where
  HVec _ x == HVec _ y = x == y

instance (Show k, H.Element k) => Show (HArr (M m n) k) where
  show (HMat _ _ x) = show x 

instance (Show k, Storable k) => Show (HArr (V n) k) where
  show (HVec _ v) = show v

instance Num (H.Matrix k) => Num (HArr (M m n) k) where
  HMat m n x + HMat m' n' y = HMat m n (x + y)
  HMat m n x - HMat m' n' y = HMat m n (x - y)
  HMat m n x * HMat m' n' y = HMat m n (x * y)
  negate (HMat m n x) = HMat m n (negate x)
  abs (HMat m n x) = HMat m n (abs x)
  signum (HMat m n x) = HMat m n (signum x)
  fromInteger = error "not allowed"


instance Num (H.Vector k) => Num (HArr (V n) k) where
  HVec n x + HVec n' y = HVec n (x + y)
  HVec n x - HVec n' y = HVec n (x - y)
  HVec n x * HVec n' y = HVec n (x * y)
  negate (HVec n x) = HVec n (negate x)
  abs (HVec n x) = HVec n (abs x)
  signum (HVec n x) = HVec n (signum x)
  fromInteger = error "not allowed"


instance H.Container H.Vector e => Scale e HArr where
  c .* HVec n v = HVec n $ H.scale c v
  c .* HMat i j m = HMat i j $ H.scale c m

instance (Floating k, Num (H.Matrix k), Num (H.Vector k), H.Field k, Numeric k) => Matr k HArr where

  toRows (HMat m n x) = MkVect (map (HVec n) (H.toRows x))
  toColumns (HMat m n x) = MkVect (map (HVec m) . H.toRows $ H.trans x)

  fromRows vect@(MkVect vs@(HVec n _ : _)) = HMat (V.length vect) n 
    (H.fromRows (map (\(HVec _ v) -> v) vs))

  fromColumns vect@(MkVect vs@(HVec m _ : _)) = HMat m (V.length vect)
    (H.fromColumns (map (\(HVec _ v) -> v) vs))

  fromVects vect@(MkVect xss@(xs : _)) = 
    HMat (V.length vect) (V.length xs) (H.fromLists (map V.toList xss))
  toVects (HMat m n x) = MkVect (map MkVect (H.toLists x))
  fromVect vect@(MkVect xs) = HVec (V.length vect) (H.fromList xs)
  toVect (HVec n v) = MkVect (H.toList v)
  takeDiag (HMat m n x) = HVec m (H.takeDiag x)

  asColMat (HVec n v) = HMat n (lit (Proxy :: Proxy 1)) (H.fromColumns [v])

  dim (HMat m n x) = (m, n)
  rows (HMat m n x) = m
  cols (HMat m n x) = n
  trans (HMat m n x) = HMat n m (H.trans x)
  ident n = HMat n n (H.ident (snat n))

  outer (HVec m u) (HVec n v) = HMat m n (H.outer u v)
  (>.<) (HVec _ u) (HVec _ v) = u H.<.> v

  mXm (HMat m n x) (HMat n' p y) = HMat m p (x H.<> y)
  chol (HMat n n' x) = HMat n' n . H.trans $ H.chol x

  fromDiag (HVec n v) = HMat n n $ H.diag v
  len (HVec n v) = n

  elementwiseprod (HMat m n x) (HMat m' n' y) = HMat m n (x * y)

  linearSolve (HMat n n' m) (HMat n'' p b) = HMat n p (H.linearSolve m b)
  trsymprod (HMat _ _ x) (HMat _ _ y) = H.sumElements (x * y)

  constant x n = HVec n (H.constant x (snat n))
