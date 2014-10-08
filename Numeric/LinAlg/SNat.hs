{-# LANGUAGE DataKinds, GADTs, KindSignatures, TypeOperators #-}
{-# LANGUAGE ScopedTypeVariables #-}
module Numeric.LinAlg.SNat (
  SNat, snat,
  lit,
  plus, times, pow,
  minus,
  cmp
) where

import Data.Proxy (Proxy (Proxy))
import Data.Type.Equality ((:~:) (Refl))
import GHC.TypeLits
import Unsafe.Coerce (unsafeCoerce)

newtype SNat :: Nat -> * where
  MkNat :: Int -> SNat n
  deriving (Eq, Show, Ord)

snat :: SNat n -> Int
snat (MkNat i) = i

lit :: KnownNat n => Proxy n -> SNat n
lit n = MkNat (fromIntegral (natVal n))

plus :: SNat n -> SNat m -> SNat (n + m)
plus (MkNat i) (MkNat j) = MkNat (i + j)

minus :: SNat n -> SNat m -> Maybe (SNat (n - m))
minus (MkNat i) (MkNat j) = let x = i - j in
  if x < 0
    then Nothing
    else Just (MkNat x)

times :: SNat n -> SNat m -> SNat (n * m)
times (MkNat i) (MkNat j) = MkNat (i * j)

pow :: SNat n -> SNat m -> SNat (n ^ m)
pow (MkNat i) (MkNat j) = MkNat (i ^ j)

cmp :: SNat n -> SNat m -> Maybe (n :~: m)
cmp (MkNat i) (MkNat j) = if i == j
  then Just (unsafeCoerce Refl)
  else Nothing

