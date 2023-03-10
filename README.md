# toms748

## C-compatible implementation of TOMS 748 root finding algorithm

This is a C-compatible implementation in C++ of TOMS 748 algorithm of Alefeld, Potra and Shi to find a root of a
function in a given interval using inverse cubic and “Newton-quadratic” interpolations.

This implementation is adapted from a combination of the SciPy Python and Boost C++ implementations.

There are two reasons for avoiding the existing Boost C++ implementation: the heavy weight of Boost which is totally
unnecessary for just using this algorithm, and that it is not C-compatible.

The C-compatibility provided here means that wrappers for this library may potentially be written in many languages.

## About the algorithm

The TOMS 748 algorithm with two interpolations per iteration is asymptotically the most efficient method known for
finding a root of a four times continuously differentiable function.

In contrast with Brent’s algorithm, which may only decrease the length of the enclosing bracket on the last step, this
algorithm decreases it at each iteration with the same asymptotic efficiency as it finds the root.

Note that this algorithm requires that the function is continuous.

Ref:
Alefeld, G. E. and Potra, F. A. and Shi, Yixun,
Algorithm 748: Enclosing Zeros of Continuous Functions,
ACM Trans. Math. Softw. Volume 221(1995) doi = {10.1145/210089.210111}

## Main Usage

To use this routine, just copy the files `toms748.cpp`, `toms748.h` and `areclose.hpp` to convenient locations in your
project, do `#include "toms748.h"` and add `toms748.cpp` in your object file compiling and linking.

The basic function is as follows:

```
double toms748(
    Toms748InputFunction function,
    void * otherInput,
    double intervalStart,
    double intervalEnd,
    Toms748ResultStatus * resultStatus);
```

`Toms748InputFunction` is just a typedef to the function pointer type `double (*)(double, void *)`.

Basically the function `function` whose root is to be searched for is passed to the root-searching routine in this
format along with the input bracket for the unknown variable in `intervalStart` and `intervalEnd`.

Additional (fixed) input data may also be provided via an opaque pointer `otherInput`, which will be passed as-is to the
function. The function is expected to cast it to a pointer to the required type and extract the additional data
therefrom.

A pointer `resultStatus` to a struct type may also be passed if more information is desired on the execution of the
routine. This information includes the number of iterations, number of function calls, actual bracket arrived at and
error code if any. Of course, this pointer may be left null if this is not desired.

The return value of the routine is the middle-point of the bracket finally arrived at, and this is nominally to be taken
as the root. This value will be `NaN` if an error was encountered.

Looking at the [C example](examples/example.c) and [C++ example](examples/example.cpp)
should make the usage much more concretely clear.

## Other Usage

The more customizable function `toms748Custom` (note the capital C) has four additional parameters at the end:

```
    double absoluteTolerance,
    double relativeTolerance,
    int maximumIterations,
    int interpolationsPerIteration
```

The `absoluteTolerance` and `relativeTolerance` are as in the [`areclose`](https://github.com/jamadagni/areclose) test
used here. The other two should be clear.

A convenience overload for `toms748` is available in C++ which can just take a functor and the input bracket. (The
`resultStatus` pointer if not provided defaults to `nullptr`.) This may make usage more conceptually convenient, but may
also mean slightly more compiled code due to the minimal template functions connecting C++ to the actual C API routines.

The convenience function `toms748d` is also available:

```
double toms748d(
    Toms748DoubleFunction function,
    double target,
    double intervalStart,
    double intervalEnd)
```

This helps search for an input value to a univariate function that will produce a given output.

Here `Toms748DoubleFunction` is just a typedef to the univariate double-valued function pointer type `double (*)(double)`.

`target` is the value that needs to be produced by `function`. `intervalStart` and `intervalEnd` are as before.

Reading the [header file](toms748.h) should also provide more insight.

## License

© 2021, Shriramana Sharma, samjnaa-at-gmail-dot-com

Use, modification and distribution are permitted subject to the "BSD-2-Clause"-type license stated in the accompanying file LICENSE.txt.
