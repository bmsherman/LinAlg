{-# LANGUAGE MultiParamTypeClasses, FunctionalDependencies #-}

{- | This module provides an embedded DSL for specifying immutable
linear algebra computations.

This is done with the 'Matr' typeclass.

-}

module Numeric.LinAlg where
import Data.List (transpose, intercalate)

-- | An instance of @'Matr' k v m@ means
-- that @v k@ is a vector type and @m k@ is a matrix type over the field @k@.
class (Mul k m m m, Mul k m v v, Mul k v m v,
       Scale k v, Scale k m,
       Solve k m v, Solve k m m,
       Num (v k), Num (m k), Floating k) 
    => Matr k v m | m -> v, v -> m where
  infixl 7 >.<

  --
  -- Data transfer
  --

  -- | Convert a list of elements to a vector.
  fromList :: [k] -> v k

  -- | Convert a vector to a list of its elements.
  toList :: v k -> [k]

  -- | Convert a row-major list of lists of elements (which should all have
  -- the same length) to a matrix containing those elements.
  fromLists :: [[k]] -> m k

  -- | Convert a matrix to a list of its rows, each given as a list of 
  -- elements.
  toLists :: m k -> [[k]]
  
  -- | Convert a matrix to a list of its rows.
  toRows :: m k -> [v k]
  toRows = map fromList . toLists

  -- | Convert a matrix to a list of its columns.
  toColumns :: m k -> [v k]
  toColumns = toRows . trans

  -- | Convert a list of vectors to a matrix having those vectors as rows.
  fromRows :: [v k] -> m k
  fromRows = fromLists . map toList

  -- | Convert a list of vectors to a matrix having those vectors as
  -- columns.

  --Is there a better default implementation for this?
  fromColumns :: [v k] -> m k
  fromColumns = trans . fromRows

  -- | Regard a vector as a matrix with a single column.
  asColMat :: v k -> m k
  asColMat = fromColumns . (:[])

  -- | Convert a matrix which has only one column to a vector.
  -- This function may have undefined behavior if the input matrix has more
  -- than one column.
  asColVec :: m k -> v k
  asColVec = head . toColumns

  -- | Produce a diagonal matrix with the given vector along its diagonal
  -- (and zeros elsewhere).
  fromDiag :: v k -> m k

  -- | Return a vector of elements along the diagonal of the matrix.
  -- Does not necessarily fail if the matrix is not square.
  takeDiag :: m k -> v k
  takeDiag = fromList . diag 0 . toLists where
    diag _ [] = []
    diag n (xs:xss) = (xs !! n) : diag (n+1) xss 

  --
  -- Core operations
  --
  
  -- | Dimension of a matrix (rows, columns).
  dim :: m k -> (Int, Int)

  -- | The number of rows in a matrix.
  rows :: m k -> Int
  rows = fst . dim

  -- | The number of columns in a matrix.
  cols :: m k -> Int
  cols = snd . dim
  
  -- | The length of a vector.
  len :: v k -> Int

  -- | Transpose a matrix.
  trans :: m k -> m k
  trans = fromLists . transpose . toLists

  -- | Construct the identity matrix of a given dimension. 
  ident :: Int -> m k
  ident n = fromLists [ [ if i==j then 1 else 0 | j <- [1..n] ] | i <- [1..n] ]

  -- | Compute the outer product of two vectors.
  outer :: v k -> v k -> m k
  
  -- | Compute the dot product (i.e., inner product) of two vectors.
  (>.<) :: v k -> v k -> k

  -- | Compute the elementwise product (i.e., Hadamard product) of 
  -- two matrices.
  elementwiseprod :: m k -> m k -> m k
  elementwiseprod x y = let [x', y'] = map toLists [x,y] in
    fromLists $ zipWith (zipWith (*)) x' y'

  --
  -- Solving linear systems
  --

  -- | General matrix inverse.
  inv :: m k -> m k
  inv l = case square l of Just n -> l <\> ident n; Nothing -> error "inv: Matrix not square."

  -- | Inverse of a lower-triangular matrix.
  invL  :: m k -> m k
  invL l = case square l of Just n -> l .\ ident n; Nothing -> error "invL: Matrix not square."

  -- | Inverse of an upper-triangular matrix.
  invU  :: m k -> m k
  invU u = case square u of Just n -> u ^\ ident n; Nothing -> error "invU: Matrix not square."


  --
  -- Functions related to Cholesky decomposition
  --

  -- | Cholesky decomposition of a positive-definite symmetric matrix.
  -- Returns lower triangular matrix of the decomposition. May not
  -- necessarily zero out the upper portion of the matrix.
  chol  :: m k -> m k

  -- | Invert a positive-definite symmetric system using a precomputed
  -- Cholesky decomposition. That is, if @ l == 'chol' a @ and 
  -- @ b == 'cholInv' l @, then @ b @ is the inverse of @ a @.
  cholInv :: m k -> m k
  cholInv l = let il = invL l in trans il >< il

  -- | Compute the log-determinant of a positive-definite symmetric matrix
  -- using its precomputed Cholesky decomposition. That is, if 
  -- @ l == 'chol' a @ and @ d == 'cholLnDet' l @, then @ d @ is the 
  -- log-determinant of @ a @.
  cholLnDet :: m k -> k
  cholLnDet = (2*) . sum . map log . toList . takeDiag

  --
  -- Other functions
  --

  -- | If matrices @ a @ and @ b @ are symmetric, then @ trsymprod a b @ is
  -- the trace of @ a '><' b @. Note that in this case, the trace of 
  -- @ a '><' b @ is just the sum of the elements in the element-wise 
  -- product of @ a @ and @ b @.
  trsymprod :: m k -> m k -> k

  -- | Create a vector of a given length whose elements all have the same
  -- value.
  constant :: k -> Int -> v k
  constant alpha n = fromList (replicate n alpha)

  


-- | If the matrix is square, return 'Just' its dimension; otherwise,
-- 'Nothing'.
square :: Matr k v m => m k -> Maybe Int
square m = let (i,j) = dim m in if i==j then Just i else Nothing

-- | Pretty-print a matrix. (Actual prettiness not guaranteed!)
showMat :: Show k => Matr k v m => m k -> String
showMat = intercalate "\n" . map ( intercalate "\t" . map show ) . toLists



-- | A typeclass representing things which can be multiplied together.
class Mul k a b c | a b -> c, a c -> b, b c -> a where
  infixl 7 ><
  (><) :: a k -> b k -> c k

-- | Scalar multiplication for vector spaces.
class Scale k v where
  infixl 7 .*
  (.*) :: k -> v k -> v k

-- | A typeclass representing things for which multiplication can possibly
-- be inverted.
class Mul k m c c => Solve k m c | c -> m where
  infixl 6 <\>
  infixl 7 ^\
  infixl 7 .\
  infixl 7 \\

  -- | Solve a general system. If there is some @x@ such that
  -- @ m '><' x == b @, then @ x == m '<\>' b @. 
  (<\>) :: m k -> c k -> c k

  -- | Solve a lower-triangular system. If @l@ is a lower-triangular 
  -- matrix, then @ x == l '.\' b @ means that @ l '><' x == b @.
  (.\) :: m k -> c k -> c k 
  (.\) = (<\>)

  -- | Solve a upper-triangular system. If @u@ is a upper-triangular
  -- matrix, then @ x == u '^\' b @ means that @ u '><' x == b @.
  (^\) :: m k -> c k -> c k
  (^\) = (<\>) 

  -- | Solve a positive-definite symmetric system. That is, if @ a @ is a
  -- positive-definite symmetric matrix and @ x == a '\\' b @, then
  -- @ a '><' x == b @.
  (\\)  :: m k -> c k -> c k
  (\\) = (<\>)
  --(\\) a = cholSolve (chol a)  
  --default implementation with Cholesky decomposition

  -- | Solve a positive-definite symmetric system using a precomputed
  -- Cholesky decomposition. If @ l == 'chol' a @ and 
  -- @ x == l \`cholSolve\` b @, then @ a '><' x == b @.
  cholSolve :: m k -> c k -> c k
  --l `cholSolve` b = trans l ^\ (l .\ b) --This may be broken? It's not working well
