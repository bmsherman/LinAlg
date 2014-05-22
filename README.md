LinAlg
======

LinAlg is a Haskell library which provides an embedded DSL for
specifying numeric linear algebra computations in a purely
functional manner. Its interface is very much in the spirit of the Haskell
library
[hmatrix](http://hackage.haskell.org/package/hmatrix).

Backends exist (packaged separately) to execute these computations on the
CPU (using [hmatrix](http://hackage.haskell.org/package/hmatrix)) or on
the GPU (using bindings to CUBLAS and MAGMA).

For example, suppose we'd like to express parameter estimation
for ordinary least-squares linear regression:

```Haskell
import Numeric.LinAlg

-- | Given a matrix of features of input data (where
-- each row is a datum) and a vector of outputs,
-- compute the parameters which minimize the sum
-- of squared error.
linearRegression :: Matr k v m => m k -> v k -> v k
linearRegression x y = (trans x >< x) <\> (trans x >< y)
```

Backends
--------

LinAlg does not come with any backends. That is, it is
impossible to *execute* any computations with this library
alone. There are two backends available:

1. **[LinAlg-hmatrix](https://github.com/bmsherman/LinAlg-hmatrix)**
for execution of computations on the CPU. It depends on 
[hmatrix](http://hackage.haskell.org/package/hmatrix), which
in turn uses a CPU BLAS library as well as LAPACK.

2. **[LinAlg-magma](https://github.com/bmsherman/LinAlg-magma)** for
execution of computations on the GPU (using CUDA). It depends
on the 
[CUBLAS](https://developer.nvidia.com/cuBLAS) and 
[MAGMA](http://icl.cs.utk.edu/magma/) libraries, and so it depends on both
[Haskell bindings to CUBLAS](https://github.com/bmsherman/cublas)
as well as
[Haskell bindings to MAGMA](https://github.com/bmsherman/magma-gpu).

### Installation

It's just a simple Cabal package:
```shell
cabal install
```
