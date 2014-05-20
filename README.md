LinAlg
======

LinAlg is a Haskell library which provides an embedded DSL for
specifying numeric linear algebra computations in a purely
functional manner.

Its interface is very much in the spirit of the Haskell library
[hmatrix](http://hackage.haskell.org/package/hmatrix).

For example, suppose we'd like to express parameter estimation
for ordinary least-squares linear regression:

```Haskell
-- | Given a matrix of features of input data (where
-- each row is a datum) and a vector of outputs,
-- compute the parameters which minimize the sum
-- of squared error.
linearRegression :: Matr m k v => m k -> v k -> v k
linearRegression x y = (trans x >< x) <\> (trans x >< y)
```

Backends
--------

LinAlg does not come with any backends. That is, it is
impossible to *execute* any computations with this library
alone. There are two backends available:

1. **LinAlg-hmatrix** for execution of computations on the CPU. It
depends on 
[hmatrix](http://hackage.haskell.org/package/hmatrix), which
in turn uses a CPU BLAS library as well as LAPACK.

2. **LinAlg-magma** for execution of computations on the GPU.
It depends on Haskell bindings to the CUBLAS and MAGMA
libraries.
