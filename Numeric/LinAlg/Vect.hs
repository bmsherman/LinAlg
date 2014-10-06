{-# LANGUAGE DataKinds, GADTs, KindSignatures, TypeOperators #-}

module Numeric.LinAlg.Vect (
  Vect (..), SomeVect (..),
  nil, cons, length,
  head,
  map, fromList, toList,
  transpose,
  generate
) where

import qualified Data.List as P
import Data.Proxy (Proxy (Proxy))

import GHC.TypeLits

import Numeric.LinAlg.SNat
import Prelude hiding (map, head, length)
import qualified Prelude as P

import Unsafe.Coerce (unsafeCoerce)

data Vect :: Nat -> * -> * where
  MkVect :: [a] -> Vect n a
  deriving (Eq, Show, Ord)

data SomeVect a where
  SomeVect :: Vect n a -> SomeVect a

length :: Vect n a -> SNat n
length (MkVect xs) = case (someNatVal (fromIntegral (P.length xs))) of
  Just (SomeNat n) -> unsafeCoerce (lit n)

map :: (a -> b) -> Vect n a -> Vect n b
map f (MkVect xs) = MkVect (P.map f xs)

nil :: Vect 0 a
nil = MkVect []

cons :: a -> Vect n a -> Vect (n + 1) a
cons x (MkVect xs) = MkVect (x : xs)

head :: (1 <= n) => Vect n a -> a
head (MkVect xs) = P.head xs

fromList :: [a] -> SomeVect a
fromList xs = SomeVect (MkVect xs) 

transpose :: Vect m (Vect n a) -> Vect n (Vect m a)
transpose (MkVect xs) = 
  MkVect (P.map MkVect (P.transpose (P.map toList xs)))

concat :: Vect m (Vect n a) -> Vect (m * n) a
concat (MkVect xs) = 
  MkVect (P.concat (P.map toList xs))

toList :: Vect n a -> [a]
toList (MkVect xs) = xs

generate :: SNat n -> (Int -> a) -> Vect n a
generate n f = MkVect (go 0) where
  n' = snat n
  go i = if i < n' then f i : go (i + 1) else []
