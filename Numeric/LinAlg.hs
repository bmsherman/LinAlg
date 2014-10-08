{-# LANGUAGE MultiParamTypeClasses, FunctionalDependencies #-}
{-# LANGUAGE DataKinds, KindSignatures, PolyKinds, TypeOperators #-}
{-# LANGUAGE GADTs #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE ScopedTypeVariables #-}

{- | This module provides an interface for specifying immutable
linear algebra computations.

This is done with the 'Matr' typeclass.
-}

module Numeric.LinAlg where

import Data.List (transpose, intercalate)
import Data.Type.Equality ((:~:))
import GHC.TypeLits (Nat)

import qualified Numeric.LinAlg.Vect as V
import Numeric.LinAlg.Vect (Vect)
import Numeric.LinAlg.SNat

data Dim = V Nat | M Nat Nat

-- | Infix, overloaded versions of functions for solving linear equations
class Solve (n :: Nat) (dim :: Dim) | dim -> n where
  infixl 7 <\>
  infixl 7 ^\
  infixl 7 .\
  infixl 7 \\

  -- | 'linearSolve'
  (<\>) :: Matr k arr => arr (M n n) k -> arr dim k -> arr dim k

  -- | 'utriSolve'
  (^\) :: Matr k arr => arr (M n n) k -> arr dim k -> arr dim k

  -- | 'ltriSolve'
  (.\) :: Matr k arr => arr (M n n) k -> arr dim k -> arr dim k

  -- | 'posdefSolve'
  (\\) :: Matr k arr => arr (M n n) k -> arr dim k -> arr dim k


instance Solve n (M n p) where
  (<\>) = linearSolve
  (^\) = utriSolve
  (.\) = ltriSolve
  (\\) = posdefSolve

instance Solve n (V n) where
  x <\> y = asColVec (linearSolve x (asColMat y))
  x ^\ y  = asColVec (utriSolve   x (asColMat y))
  x .\ y  = asColVec (ltriSolve   x (asColMat y))
  x \\ y  = asColVec (posdefSolve x (asColMat y))

class Product (a :: Dim) (b :: Dim) where
  type Prod a b :: Dim
  infixl 7 ><
  (><) :: Matr k arr => arr a k -> arr b k -> arr (Prod a b) k

instance Product (M m n) (M n p) where
  type Prod (M m n) (M n p) = M m p
  x ><  y = mXm x y

instance Product (V  n) (M n p) where
  type Prod (V   n) (M n p) = V   p
  x >< y = asColVec (trans (mXm (trans (asColMat x)) y))

instance Product (M m n) (V n  ) where
  type Prod (M m n) (V n  ) = V m 
  x >< y = asColVec (mXm x (asColMat y))

-- | An instance of @'Matr' k v m@ means
-- that @v k@ is a vector type and @m k@ is a matrix type over the field @k@.
class  ( Floating k, Scale k arr) 
  => Matr (k :: *) (arr :: Dim -> * -> *) where
  infixl 7 >.<

  --
  -- Data transfer
  --

  -- | Convert a list of elements to a vector.
  fromVect :: Vect n k -> arr (V n) k

  -- | Convert a vector to a list of its elements.
  toVect :: arr (V n) k -> Vect n k

  -- | Convert a row-major list of lists of elements (which should all have
  -- the same length) to a matrix containing those elements.
  fromVects :: Vect m (Vect n k) -> arr (M m n) k

  -- | Convert a matrix to a list of its rows, each given as a list of 
  -- elements.
  toVects :: arr (M m n) k -> Vect m (Vect n k)
  
  -- | Convert a matrix to a list of its rows.
  toRows :: arr (M m n) k -> Vect m (arr (V n) k)

  -- | Convert a matrix to a list of its columns.
  toColumns :: arr (M m n) k -> Vect n (arr (V m) k)
  toColumns = toRows . trans

  -- | Convert a list of vectors to a matrix having those vectors as rows.
  fromRows :: Vect m (arr (V n) k) -> arr (M m n) k

  -- | Convert a list of vectors to a matrix having those vectors as
  -- columns.
  fromColumns :: Vect n (arr (V m) k) -> arr (M m n) k

  -- | Regard a vector as a matrix with a single column.
  asColMat :: arr (V n) k -> arr (M n 1) k

  -- | Convert a matrix which has only one column to a vector.
  -- This function may have undefined behavior if the input matrix has more
  -- than one column.
  asColVec :: arr (M n 1)  k -> arr (V n) k
  asColVec = V.head . toColumns

  -- | Produce a diagonal matrix with the given vector along its diagonal
  -- (and zeros elsewhere).
  fromDiag :: arr (V n) k -> arr (M n n) k

  -- | Return a vector of elements along the diagonal of the matrix.
  -- Does not necessarily fail if the matrix is not square.
  takeDiag :: arr (M n n) k -> arr (V n) k

  --
  -- Core operations
  --
  
  -- | Dimension of a matrix (rows, columns).
  dim :: arr (M m n) k -> (SNat m, SNat n)

  -- | The number of rows in a matrix.
  rows :: arr (M m n) k -> SNat m
  rows = fst . dim

  -- | The number of columns in a matrix.
  cols :: arr (M m n) k -> SNat n
  cols = snd . dim
  
  -- | The length of a vector.
  len :: arr (V n) k -> SNat n

  -- | Transpose a matrix.
  trans :: arr (M m n) k -> arr (M n m) k

  -- | Construct the identity matrix of a given dimension. 
  ident :: SNat n -> arr (M n n) k

  -- | Compute the outer product of two vectors.
  outer :: arr (V m) k -> arr (V n) k -> arr (M m n) k
  
  -- | Compute the dot product (i.e., inner product) of two vectors.
  (>.<) :: arr (V n) k -> arr (V n) k -> k

  -- Multiplication
  mXm :: arr (M m n) k -> arr (M n p) k -> arr (M m p) k

  -- | Compute the elementwise product (i.e., Hadamard product) of 
  -- two matrices.
  elementwiseprod :: arr (M m n) k -> arr (M m n) k -> arr (M m n) k

  --
  -- Solving linear systems
  --

  -- | General matrix inverse.
  inv :: arr (M n n) k -> arr (M n n) k
  inv m = linearSolve m (ident (rows m))

  -- | Inverse of a lower-triangular matrix.
  invL  :: arr (M n n) k -> arr (M n n) k
  invL = inv

  -- | Inverse of an upper-triangular matrix.
  invU  :: arr (M n n) k -> arr (M n n) k
  invU = inv


  --
  -- Functions related to Cholesky decomposition
  --

  -- | Cholesky decomposition of a positive-definite symmetric matrix.
  -- Returns lower triangular matrix of the decomposition. May not
  -- necessarily zero out the upper portion of the matrix.
  chol  :: arr (M n n) k -> arr (M n n) k

  -- | Invert a positive-definite symmetric system using a precomputed
  -- Cholesky decomposition. That is, if @ l == 'chol' a @ and 
  -- @ b == 'cholInv' l @, then @ b @ is the inverse of @ a @.
  cholInv :: arr (M n n) k -> arr (M n n) k
  cholInv l = let il = invL l in trans il >< il

  -- | Compute the log-determinant of a positive-definite symmetric matrix
  -- using its precomputed Cholesky decomposition. That is, if 
  -- @ l == 'chol' a @ and @ d == 'cholLnDet' l @, then @ d @ is the 
  -- log-determinant of @ a @.
  cholLnDet :: arr (M n n) k -> k
  cholLnDet = (2*) . sum . map log . V.toList . toVect . takeDiag

  --
  -- Other functions
  --

  -- | If matrices @ a @ and @ b @ are symmetric, then @ trsymprod a b @ is
  -- the trace of @ a '><' b @. Note that in this case, the trace of 
  -- @ a '><' b @ is just the sum of the elements in the element-wise 
  -- product of @ a @ and @ b @.
  trsymprod :: arr (M n n) k -> arr (M n n) k -> k

  -- | Create a vector of a given length whose elements all have the same
  -- value.
  constant :: k -> SNat n -> arr (V n) k

  --
  -- Solving linear equations
  --

  -- | Solve a general system. If there is some @x@ such that
  -- @ m '><' x == b @, then @ x == 'linearSolve m b' @. 
  linearSolve :: arr (M n n) k -> arr (M n p) k -> arr (M n p) k

  -- | Solve a lower-triangular system. If @l@ is a lower-triangular 
  -- matrix, then @ x == 'ltriSolve' l b @ means that @ l '><' x == b @.
  ltriSolve :: arr (M n n) k -> arr (M n p) k -> arr (M n p) k
  ltriSolve = linearSolve

  -- | Solve a upper-triangular system. If @u@ is a upper-triangular
  -- matrix, then @ x == 'utriSolve' u b @ means that @ u '><' x == b @.
  utriSolve :: arr (M n n) k -> arr (M n p) k -> arr (M n p) k
  utriSolve = linearSolve

  -- | Solve a positive-definite symmetric system. That is, if @ a @ is a
  -- positive-definite symmetric matrix and @ x == 'posdefSolve' a b @, then
  -- @ a '><' x == b @.
  posdefSolve  :: arr (M n n) k -> arr (M n p) k -> arr (M n p) k
  posdefSolve = linearSolve

  -- | Solve a positive-definite symmetric system using a precomputed
  -- Cholesky decomposition. If @ l == 'chol' a @ and 
  -- @ x == l \`'cholSolve'\` b @, then @ a '><' x == b @.
  cholSolve :: arr (M n n) k -> arr (M n p) k -> arr (M n p) k
  l `cholSolve` b = trans l ^\ (l .\ b) --This may be broken? It's not working well

-- | If the matrix is square, return 'Just' its dimension; otherwise,
-- 'Nothing'.
square :: Matr k arr => arr (M m n) k -> Maybe (m :~: n)
square m = let (i,j) = dim m in cmp i j

-- | Pretty-print a matrix. (Actual prettiness not guaranteed!)
showMat :: (Show k, Matr k arr) => arr (M m n) k -> String
showMat = intercalate "\n" . map ( intercalate "\t" . map show ) . toLists


toLists :: Matr k arr => arr (M m n) k -> [[k]]
toLists = map (V.toList) . V.toList . toVects

-- | Scalar multiplication for vector spaces.
class Scale k arr where
  infixl 7 .*
  (.*) :: k -> arr dim k -> arr dim k
