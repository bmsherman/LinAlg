LinAlg
======

LinAlg is a Haskell library which provides an interface for
specifying numeric linear algebra computations in a purely
functional manner. Its interface is very much in the spirit of the Haskell
library
[hmatrix](http://hackage.haskell.org/package/hmatrix).
Array sizes are reflected statically in the types.

Backends exist to execute these computations on the
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
linearRegression :: Matr k arr => arr (M n m) k -> arr (V n) k -> arr (V m) k
linearRegression x y = (trans x >< x) <\> (trans x >< y)
```

Backends
--------

1) **HMatrix** (included in this package)
for execution of computations on the CPU. It depends on 
[hmatrix](http://hackage.haskell.org/package/hmatrix), which
in turn uses a CPU BLAS library as well as LAPACK.

It's a rather straightforward translation from hmatrix to LinAlg; in fact,
the names for many functions are exactly the same.

The HMatrix backend is installed by default, but can be disabled
by disabling the "hmatrix" flag (so the dependency on the hmatrix package
is obviated).

2) **[LinAlg-magma](https://github.com/bmsherman/LinAlg-magma)** 
(an external package)
for execution of computations on the GPU (using CUDA). It depends
on the 
[CUBLAS](https://developer.nvidia.com/cuBLAS) and 
[MAGMA](http://icl.cs.utk.edu/magma/) libraries, and so it depends on both
[Haskell bindings to CUBLAS](https://github.com/bmsherman/cublas)
as well as
[Haskell bindings to MAGMA](https://github.com/bmsherman/magma-gpu).

### Documentation

[See the Haddock documentation](http://bmsherman.github.io/haddock/LinAlg/index.html).

### Installation

It's just a simple Cabal package:
```shell
cabal install
```
