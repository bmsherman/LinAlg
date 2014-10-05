{-# LANGUAGE MultiParamTypeClasses, FlexibleInstances, UndecidableInstances #-}
{-# LANGUAGE DataKinds, KindSignatures, GADTs #-}
{-# LANGUAGE ScopedTypeVariables #-}

{- | A straightforward LinAlg backend using the hmatrix library. -}

module Numeric.LinAlg.HMatrix (
    module Numeric.LinAlg
  , HArr (HMat, HVec)
) where

import Numeric.LinAlg
import qualified Numeric.LinearAlgebra as H
import Numeric.LinearAlgebra.HMatrix (Numeric)

import Foreign.Storable (Storable)

data HArr :: Dim -> * -> * where
  HMat :: H.Matrix e -> HArr (M m n) e
  HVec :: H.Vector e -> HArr (V n  ) e


instance (H.Container H.Vector k) => Eq (HArr (M m n) k) where
  HMat x == HMat y = x == y


instance (Eq k, Storable k) => Eq (HArr (V n) k) where
  HVec x == HVec y = x == y

instance (Show k, H.Element k) => Show (HArr (M m n) k) where
  show (HMat x) = show x 

instance (Show k, Storable k) => Show (HArr (V n) k) where
  show (HVec v) = show v

instance Num (H.Matrix k) => Num (HArr (M m n) k) where
  HMat x + HMat y = HMat (x + y)
  HMat x - HMat y = HMat (x - y)
  HMat x * HMat y = HMat (x * y)
  negate (HMat x) = HMat (negate x)
  abs (HMat x) = HMat (abs x)
  signum (HMat x) = HMat (signum x)
  fromInteger = HMat . fromInteger


instance Num (H.Vector k) => Num (HArr (V n) k) where
  HVec x + HVec y = HVec (x + y)
  HVec x - HVec y = HVec (x - y)
  HVec x * HVec y = HVec (x * y)
  negate (HVec x) = HVec (negate x)
  abs (HVec x) = HVec (abs x)
  signum (HVec x) = HVec (signum x)
  fromInteger = HVec . fromInteger


instance H.Container H.Vector e => Scale e HArr where
  c .* HVec v = HVec $ H.scale c v
  c .* HMat m = HMat $ H.scale c m

instance (Floating k, Num (H.Matrix k), Num (H.Vector k), H.Field k, Numeric k) => Matr k HArr where

  toRows (HMat m) = map HVec (H.toRows m)
  toColumns (HMat m) = map HVec . H.toRows $ H.trans m
  fromRows = HMat . H.fromRows . map (\(HVec v) -> v)
  fromColumns = HMat . H.fromColumns . map (\(HVec v) -> v)

  fromLists = HMat . H.fromLists
  toLists (HMat m) = H.toLists m
  fromList = HVec . H.fromList
  toList (HVec v) = H.toList v
  takeDiag (HMat m) = HVec (H.takeDiag m)

  dim (HMat m) = (H.rows m, H.cols m)
  rows (HMat m) = H.rows m
  cols (HMat m) = H.cols m
  trans (HMat m) = HMat (H.trans m)
  ident = HMat . H.ident

  outer (HVec u) (HVec v) = HMat (H.outer u v)
  (>.<) (HVec u) (HVec v) = u H.<.> v

  mXm (HMat x) (HMat y) = HMat (x H.<> y)
  chol (HMat x) = HMat . H.trans $ H.chol x

  fromDiag (HVec v) = HMat $ H.diag v
  len (HVec v) = H.dim v

  elementwiseprod (HMat x) (HMat y) = HMat (x * y)

  linearSolve (HMat m) (HMat b) = HMat (H.linearSolve m b)
  trsymprod (HMat x) (HMat y) = H.sumElements (x * y)

  constant x = HVec . H.constant x
