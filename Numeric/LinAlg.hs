{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE DataKinds, KindSignatures, PolyKinds #-}
{-# LANGUAGE TypeFamilies #-}

{- | This module provides an embedded DSL for specifying immutable
linear algebra computations.

This is done with the 'Matr' typeclass.

-}

module Numeric.LinAlg where
import Data.List (transpose, intercalate)
import GHC.TypeLits (Nat)

data Dim = V Nat | M Nat Nat

class Product (a :: Dim) (b :: Dim) where
  type Prod a b :: Dim
  (><) :: Matr k arr => arr a k -> arr b k -> arr (Prod a b) k

instance Product (M m n) (M n p) where
  type Prod (M m n) (M n p) = M m p
  x >< y = mXm x y

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
  -- infixl 7 >.<

  --
  -- Data transfer
  --

  -- | Convert a list of elements to a vector.
  fromList :: [k] -> arr (V n) k

  -- | Convert a vector to a list of its elements.
  toList :: arr (V n) k -> [k]

  -- | Convert a row-major list of lists of elements (which should all have
  -- the same length) to a matrix containing those elements.
  fromLists :: [[k]] -> arr (M m n) k

  -- | Convert a matrix to a list of its rows, each given as a list of 
  -- elements.
  toLists :: arr (M m n) k -> [[k]]
  
  -- | Convert a matrix to a list of its rows.
  toRows :: arr (M m n) k -> [arr (V n) k]
  toRows = map fromList . toLists

  -- | Convert a matrix to a list of its columns.
  toColumns :: arr (M m n) k -> [arr (V m) k]
  toColumns = toRows . trans

  -- | Convert a list of vectors to a matrix having those vectors as rows.
  fromRows :: [arr (V n) k] -> arr (M m n) k
  fromRows = fromLists . map toList

  -- | Convert a list of vectors to a matrix having those vectors as
  -- columns.

  --Is there a better default implementation for this?
  fromColumns :: [arr (V m) k] -> arr (M m n) k
  fromColumns = trans . fromRows

  -- | Regard a vector as a matrix with a single column.
  asColMat :: arr (V n) k -> arr (M n 1) k
  asColMat = fromColumns . (:[])

  -- | Convert a matrix which has only one column to a vector.
  -- This function may have undefined behavior if the input matrix has more
  -- than one column.
  asColVec :: arr (M n 1)  k -> arr (V n) k
  asColVec = head . toColumns

  -- | Produce a diagonal matrix with the given vector along its diagonal
  -- (and zeros elsewhere).
  fromDiag :: arr (V n) k -> arr (M n n) k

  -- | Return a vector of elements along the diagonal of the matrix.
  -- Does not necessarily fail if the matrix is not square.
  takeDiag :: arr (M n n) k -> arr (V n) k
  takeDiag = fromList . diag 0 . toLists where
    diag _ [] = []
    diag n (xs:xss) = (xs !! n) : diag (n+1) xss 

  --
  -- Core operations
  --
  
  -- | Dimension of a matrix (rows, columns).
  dim :: arr (M m n) k -> (Int, Int)

  -- | The number of rows in a matrix.
  rows :: arr (M m n) k -> Int
  rows = fst . dim

  -- | The number of columns in a matrix.
  cols :: arr (M m n) k -> Int
  cols = snd . dim
  
  -- | The length of a vector.
  len :: arr (V n) k -> Int

  -- | Transpose a matrix.
  trans :: arr (M m n) k -> arr (M n m) k
  trans = fromLists . transpose . toLists

  -- | Construct the identity matrix of a given dimension. 
  ident :: Int -> arr (M n n) k
  ident n = fromLists [ [ if i==j then 1 else 0 | j <- [1..n] ] | i <- [1..n] ]

  -- | Compute the outer product of two vectors.
  outer :: arr (V m) k -> arr (V n) k -> arr (M m n) k
  
  -- | Compute the dot product (i.e., inner product) of two vectors.
  (>.<) :: arr (V n) k -> arr (V n) k -> k

  -- Multiplication
  mXm :: arr (M m n) k -> arr (M n p) k -> arr (M m p) k

  -- | Compute the elementwise product (i.e., Hadamard product) of 
  -- two matrices.
  elementwiseprod :: arr (M m n) k -> arr (M m n) k -> arr (M m n) k
  elementwiseprod x y = let [x', y'] = map toLists [x,y] in
    fromLists $ zipWith (zipWith (*)) x' y'

  --
  -- Solving linear systems
  --

  -- | General matrix inverse.
  inv :: arr (M n n) k -> arr (M n n) k
  inv l = case square l of Just n -> l <\> ident n; Nothing -> error "inv: Matrix not square."

  -- | Inverse of a lower-triangular matrix.
  invL  :: arr (M n n) k -> arr (M n n) k
  invL l = case square l of Just n -> l .\ ident n; Nothing -> error "invL: Matrix not square."

  -- | Inverse of an upper-triangular matrix.
  invU  :: arr (M n n) k -> arr (M n n) k
  invU u = case square u of Just n -> u ^\ ident n; Nothing -> error "invU: Matrix not square."


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
  cholLnDet = (2*) . sum . map log . toList . takeDiag

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
  constant :: k -> Int -> arr (V n) k
  constant alpha n = fromList (replicate n alpha)

  --
  -- Solving linear equations
  --
  infixl 6 <\>
  infixl 7 ^\
  infixl 7 .\
  infixl 7 \\

  -- | Solve a general system. If there is some @x@ such that
  -- @ m '><' x == b @, then @ x == m '<\>' b @. 
  (<\>) :: arr (M m n) k -> arr (M m p) k -> arr (M m p) k

  -- | Solve a lower-triangular system. If @l@ is a lower-triangular 
  -- matrix, then @ x == l '.\' b @ means that @ l '><' x == b @.
  (.\) :: arr (M m n) k -> arr (M m p) k -> arr (M m p) k
  (.\) = (<\>)

  -- | Solve a upper-triangular system. If @u@ is a upper-triangular
  -- matrix, then @ x == u '^\' b @ means that @ u '><' x == b @.
  (^\) :: arr (M m n) k -> arr (M m p) k -> arr (M m p) k
  (^\) = (<\>) 

  -- | Solve a positive-definite symmetric system. That is, if @ a @ is a
  -- positive-definite symmetric matrix and @ x == a '\\' b @, then
  -- @ a '><' x == b @.
  (\\)  :: arr (M m n) k -> arr (M m p) k -> arr (M m p) k
  (\\) = (<\>)
  --(\\) a = cholSolve (chol a)  
  --default implementation with Cholesky decomposition

  -- | Solve a positive-definite symmetric system using a precomputed
  -- Cholesky decomposition. If @ l == 'chol' a @ and 
  -- @ x == l \`cholSolve\` b @, then @ a '><' x == b @.
  cholSolve :: arr (M m n) k -> arr (M m p) k -> arr (M m p) k
  --l `cholSolve` b = trans l ^\ (l .\ b) --This may be broken? It's not working well

-- | If the matrix is square, return 'Just' its dimension; otherwise,
-- 'Nothing'.
square :: Matr k arr => arr (M m n) k -> Maybe Int
square m = let (i,j) = dim m in if i==j then Just i else Nothing

-- | Pretty-print a matrix. (Actual prettiness not guaranteed!)
showMat :: (Show k, Matr k arr) => arr (M m n) k -> String
showMat = intercalate "\n" . map ( intercalate "\t" . map show ) . toLists



-- | Scalar multiplication for vector spaces.
class Scale k arr where
  infixl 7 .*
  (.*) :: k -> arr dim k -> arr dim k
