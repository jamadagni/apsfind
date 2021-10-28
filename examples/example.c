#include "../toms748.h"
#include "elliptic_integral.h"
#include <math.h>
#include <stdio.h>

// adapted from https://www.boost.org/doc/libs/1_77_0/libs/math/doc/html/math_toolkit/root_finding_examples/elliptic_eg.html
// elliptic_integral code from https://people.sc.fsu.edu/~jburkardt/c_src/elliptic_integral/elliptic_integral.html

typedef struct { double arc, rad; } _RootFunctionParams;

double _rootFunction(double unknownRadius, void * otherInput)
{
    _RootFunctionParams * p = otherInput;
    double a = fmax(unknownRadius, p->rad);
    double b = fmin(unknownRadius, p->rad);
    double k = sqrt(1 - (b * b) / (a * a));
    return 4 * a * elliptic_ek(k) - p->arc;
}

double ellipticRoot(double arcLength, double knownRadius)
{
    double guess = sqrt(arcLength * arcLength / 16 - knownRadius * knownRadius);
    _RootFunctionParams params;
    params.arc = arcLength;
    params.rad = knownRadius;
    return toms748(_rootFunction, (void *)(&params), guess - arcLength / 10, guess + arcLength / 10, NULL);
}

int main()
{
    printf("%10.6f\n", ellipticRoot(300, 28));
}
