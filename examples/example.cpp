#include "../toms748.h"
#include "elliptic_integral.h"
#include <cmath>
#include <cstdio>

// adapted from https://www.boost.org/doc/libs/1_77_0/libs/math/doc/html/math_toolkit/root_finding_examples/elliptic_eg.html
// elliptic_integral code from https://people.sc.fsu.edu/~jburkardt/c_src/elliptic_integral/elliptic_integral.html

struct RootFunctor
{
    RootFunctor(double arcLength, double knownRadius): arc{arcLength}, rad{knownRadius} {}
    double arc, rad;

    double operator()(double unknownRadius)
    {
        double a = fmax(unknownRadius, rad);
        double b = fmin(unknownRadius, rad);
        double k = sqrt(1 - (b * b) / (a * a));
        return 4 * a * elliptic_ek(k) - arc;
    }
};

double ellipticRoot(double arcLength, double knownRadius)
{
    double guess = sqrt(arcLength * arcLength / 16 - knownRadius * knownRadius);
    RootFunctor functor{arcLength, knownRadius};
    return toms748(functor, guess - arcLength / 10, guess + arcLength / 10, nullptr);
}

int main()
{
    printf("%10.6f\n", ellipticRoot(300, 28));
}
