# apsfind

### C-compatible implementation of Alefeld, Potra and Shi's root finding algorithm (TOMS 748)

## About the algorithm

Alefeld, Potra and Shi devised an algorithm to find a root of a function in a given interval using inverse cubic and
“Newton-quadratic” interpolations. It was published as Algorithm 748 of the “Transactions on Mathematical Software”,
and is hence sometimes referred to as TOMS 748.

This algorithm with two interpolations per iteration is asymptotically the most efficient method known for
finding a root of a four times continuously differentiable function.

In contrast with Brent’s algorithm, which may only decrease the length of the enclosing bracket on the last step, this
algorithm decreases it at each iteration with the same asymptotic efficiency as it finds the root.

Note that this algorithm requires that the function is continuous.

**Ref:**
Alefeld, G. E. and Potra, F. A. and Shi, Yixun.
*Algorithm 748: Enclosing Zeros of Continuous Functions.*
ACM Trans. Math. Softw. Volume 221(1995) doi = {10.1145/210089.210111}

## Present implementation

This is a C-compatible implementation of this algorithm in C++. It is adapted from a combination of the [SciPy Python](https://github.com/scipy/scipy/blob/main/scipy/optimize/_zeros_py.py) and [Boost C++](https://github.com/boostorg/math/blob/master/include/boost/math/tools/toms748_solve.hpp) implementations.

There are two reasons for avoiding the existing Boost C++ implementation: the heavy weight of Boost which is totally
unnecessary for just using this algorithm, and that it is not C-compatible.

The C-compatibility provided here means that wrappers for this library may potentially be written in many languages.

The library and its functions were earlier named `toms748` like in the source SciPy and Boost implementations. However now they are named `apsfind` to honour the authors of the algorithm, as in “Newton's method”, and for readability, as in “who is Tom and what's that 748 in there?”.

## Main Usage

To use this routine, just copy the files `apsfind.cpp`, `apsfind.h` and `areclose.hpp` to convenient locations in your
project, do `#include "apsfind.h"` and add `apsfind.cpp` in your object file compiling and linking.

The basic function is as follows:

```
double apsfind(
    ApsFindInputFunction fn,
    void * otherInput,
    double intervalStart,
    double intervalEnd,
    ApsFindResultStatus * stat)
```

`ApsFindInputFunction` is just a typedef to the function pointer type `double (*)(double, void *)`.

Basically the function `fn` whose root is to be searched for is passed to the root-searching routine in this
format along with the input interval for the unknown variable as `intervalStart` and `intervalEnd`.

Additional (known) input data may also be provided via an opaque pointer `otherInput`, which will be passed as-is to the
function. The function is expected to cast it to a pointer to the required type and extract the additional data
therefrom.

A pointer `stat` to a struct type may also be passed if more information is desired on the execution of the
routine. This information includes the number of iterations, number of function calls, actual bracket arrived at and
error code if any. Of course, this pointer may be left null if this is not desired.

The return value of the routine is the middle-point of the bracket finally arrived at, and this is nominally to be taken
as the root. This value will be `NaN` if an error was encountered.

Looking at the [C example](examples/example.c) and [C++ example](examples/example.cpp)
should make the usage much more concretely clear.

## Other Usage

### `apsfindCustom`

This more customizable function (note the capital C in the name) has four additional parameters at the end:

```
double apsfindCustom(
…
    double absoluteTolerance,
    double relativeTolerance,
    int maximumIterations,
    int interpolationsPerIteration)
```

The `absoluteTolerance` and `relativeTolerance` are as in the [`areclose`](https://github.com/jamadagni/areclose) test
used here. The other two should be clear.

### `apsfindu`

This helps search for an input value to a univariate function (hence the “u”) that will produce a given output:

```
double apsfindu(
    ApsFindUniFunction fn,
    double target,
    double intervalStart,
    double intervalEnd)
```

`ApsFindUniFunction` is just a typedef to the univariate double-valued function pointer type `double (*)(double)`.

`target` is the value that needs to be produced by `fn`. The other two are as before.

### `apsfindResultStatusPrint`

This function can help with a detailed description of the result status:

```
void apsfindResultStatusPrint(FILE * file, ApsFindResultStatus stat, int precision)
```

The `precision` is the number of significant digits used to display the bracket.

### C++ overload of `apsfind`

A convenience overload for `apsfind` is available in C++ which can just take a functor and the input interval:

```
double apsfind(
    Functor fr,
    double intervalStart,
    double intervalEnd,
    ApsFindResultStatus * stat = nullptr)
```

### Other

Reading the [header file](apsfind.h) should also provide more insight.

Also note that the prefix `apsfind` to all additional functions is used *without* CamelCase for consistency with the simple
all-small form of the main function.

## License

© 2021, Shriramana Sharma, samjnaa-at-gmail-dot-com

Use, modification and distribution are permitted subject to the "BSD-2-Clause"-type license stated in the accompanying file LICENSE.txt.
