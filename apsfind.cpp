// apsfind
// =======
//
// The algorithm of Alefeld, Potra and Shi to find a root of a function in
// a given interval using inverse cubic and “Newton-quadratic” interpolations.
// This was published as TOMS 748; see reference below.
//
// This algorithm with two interpolations per iteration is asymptotically the
// most efficient method known for finding a root of a four times continuously
// differentiable function.
//
// In contrast with Brent’s algorithm, which may only decrease the length
// of the enclosing bracket on the last step, this algorithm decreases it
// at each iteration with the same asymptotic efficiency as it finds the root.
//
// NOTE: This algorithm requires that the function is continuous.
//
// Ref:
// Alefeld, G. E. and Potra, F. A. and Shi, Yixun,
// Algorithm 748: Enclosing Zeros of Continuous Functions,
// ACM Trans. Math. Softw. Volume 221(1995) doi = {10.1145/210089.210111}
//
// *************************************************************************************************
//
// This is a C++ implementation with a C-compatible interface and an optional C++ helper.
// Copyright © 2021, Shriramana Sharma, samjnaa-at-gmail-dot-com
//
// Use, modification and distribution are permitted subject to the
// "BSD-2-Clause"-type license stated in the accompanying file LICENSE.txt
//
// *************************************************************************************************
//
// This implementation is adapted from a combination of the SciPy and Boost implementations:
// SciPy: https://github.com/scipy/scipy/blob/main/scipy/optimize/_zeros_py.py
// Boost: https://github.com/boostorg/math/blob/master/include/boost/math/tools/toms748_solve.hpp
// Parts not marked particularly either way are common to both SciPy and Boost.
//
// The SciPy implementation is licensed permissively under a “BSD-3-Clause” license.
// The Boost implementation is licensed permissively under the Boost license.
// This implementation is likewise licensed permissively under a “BSD-2-Clause” license.
//
// *************************************************************************************************

#include "apsfind.h"
#include "areclose.hpp"
#include <cmath>
#include <cfloat>
using namespace std;


// default values below from SciPy, except as marked
static const int
    defaultInterpolationsPerIteration = 2,  // from Boost, since it's the most efficient
    defaultMaximumIterations = 100;
static const double
    eps = DBL_EPSILON,
    mu = 0.5,  // factor by which brackets are expected to decrease per iteration
    defaultRelativeTolerance = eps * 4,
    defaultAbsoluteTolerance = 2e-12;


// The idea of the generic pointer for other fixed input is from:
// SciPy Cython: https://docs.scipy.org/doc/scipy/reference/optimize.cython_optimize.html
//
// The ResultStatus struct combines SciPy's RootResults class and Boost's returning the final bracket
//
// Just as a matter of internal convention, “interval” is used only for input, and when one has first
// hand checked that it “bracket”-s a root, it is further referred to as a “bracket”.

struct ApsFindSolver
{
private:
    ApsFindInputFunction f;
    void * args;
    double a, fa, b, fb, c, fc, d, fd, e, fe, u, fu;
    ApsFindResultStatus * stat;
    AreClose areClose;
    int maxiter, remiter, fncount, k;
    double root;

public:
    ApsFindSolver(
        ApsFindInputFunction fn,
        void * otherInput,
        double intervalStart,
        double intervalEnd,
        ApsFindResultStatus * resultStatus,
        double absoluteTolerance,
        double relativeTolerance,
        int maximumIterations,
        int interpolationsPerIteration):
        f{fn},
        args{otherInput},
        a{intervalStart},
        b{intervalEnd},
        stat{resultStatus},
        areClose{absoluteTolerance, relativeTolerance},
        maxiter{maximumIterations},
        remiter{maximumIterations},
        k{interpolationsPerIteration},
        fncount{0},
        e{NAN},  // to ensure first interpolation of second iteration is quadratic
        fd{0},  // dummy value to fix weird error with optimization flag
        fe{0}
    {
        if (checkBracketForEnd())  // maybe already input is so
            return;
        if (callFunctionAndCheckEnd(a, fa))
            return;
        if (callFunctionAndCheckEnd(b, fb))
            return;
        if (dontBracketARoot(fa, fb))
        {
            finalize(NAN, APSFIND_INTERVAL_DOES_NOT_BRACKET_A_ROOT, a, b);
            return;
        }

        // first iteration

        interpolateSecant();
        if (evalRebracketAndCheckEnd())
            return;
        --remiter;

        // iteration loop

        for (; remiter; --remiter)
        {
            double initBracketWidth = b - a;
            for (int stepCount = 2; stepCount < k + 2; ++stepCount)
            {
                if (cantDoCubic())
                    interpolateNewtonQuadratic(stepCount);
                else
                {
                    interpolateCubic();
                    if (outOfBracket(c))
                        interpolateNewtonQuadratic(stepCount);
                }
                if (evalRebracketAndCheckEnd())
                    return;
            }

            interpolateDoubleLengthSecant();
            if (evalRebracketAndCheckEnd())
                return;

            // bisect if bracket width did not decrease enough
            if (b - a > mu * initBracketWidth)
                interpolateBisection();
            if (evalRebracketAndCheckEnd())
                return;
        }

        finalize(NAN, APSFIND_MAXIMUM_ITERATIONS_REACHED, a, b);
    }

    double result() const { return root; }

private:
    void finalize(double val, int errorCode, double start, double end)
    {
        root = val;
        if (stat)
        {
            stat->iterations = maxiter - remiter;
            stat->functionCalls = fncount;
            stat->errorCode = errorCode;
            stat->bracketStart = start;
            stat->bracketEnd = end;
        }
    }

    bool checkBracketForEnd()
    {
        if (!areClose(a, b))
            return false;

        // final root value returned as midpoint of bracket
        interpolateBisection();
        finalize(c, APSFIND_NO_ERROR, a, b);
        return true;
    }

    bool callFunctionAndCheckEnd(double x, double & fx)
    {
        fx = f(x, args);
        fncount += 1;
        if (fx == 0)
        {
            finalize(x, APSFIND_NO_ERROR, x, x);
            return true;
        }
        if (!isfinite(fx))
        {
            finalize(NAN, APSFIND_INVALID_FUNCTION_VALUE, a, b);
            return true;
        }
        return false;
    }

    static bool dontBracketARoot(double p, double q)
    {
        return p * q > 0;
    }

    bool outOfBracket(double val)
    {
        return val <= a || val >= b;
    }

    void interpolateBisection()
    {
        // from Boost; SciPy just uses (a + b) / 2
        c = a + (b - a) / 2;
    }

    void interpolateSecant()
    {
        // from SciPy
#if 0
        // This test does not make sense because interpolation is only ever called
        // if the interval brackets a root
        if (fa == fb)
            c = NAN;
        else
#endif
        if (fabs(fb) > fabs(fa))
            c = (a - fa / fb * b) / (1 - fa / fb);
        else
            c = (b - fb / fa * a) / (1 - fb / fa);

        // eps-based test below from Boost also catches when c sits on the edges of
        // or falls outside the bracket; SciPy only has ordinary outOfBracket test
        // everywhere before doing bisection. Boost also uses and justifies this only
        // for secant but not for double-length secant or elsewhere. Not sure why.
        if (c < a + fabs(a) * eps * 5 ||
            c > b - fabs(b) * eps * 5)
            interpolateBisection();
    }

    void interpolateDoubleLengthSecant()
    {
        if (fabs(fa) < fabs(fb))
        {
            u = a; fu = fa;
        }
        else
        {
            u = b; fu = fb;
        }
        c = u - 2 * (fu / (fb - fa)) * (b - a);

        double sepCU = fabs(c - u);
        if (sepCU > (b - a) / 2)
            interpolateBisection();
        else if (sepCU < eps * u)  // from SciPy
        {
            int aExp, bExp, largeExp, smallExp;
            double nonU;
            frexp(a, & aExp);
            frexp(b, & bExp);
            if (aExp < bExp)
            {
                smallExp = aExp; largeExp = bExp; nonU = b;
            }
            else
            {
                smallExp = bExp; largeExp = aExp; nonU = a;
            }
            if (smallExp < largeExp - 50)
                c = (31 * u + nonU) / 32;
            else
            {
                double adj = fmax(areClose.absoluteTolerance(), areClose.relativeTolerance() * fabs(c));
                if (u == a)
                    c = u + adj;
                else
                    c = u - adj;
                if (outOfBracket(c))
                    interpolateBisection();
            }
        }
    }

    static double safeDiv(double n, double d, double r)
    {
        // from Boost
        if (fabs(d) < 1 && fabs(d) * DBL_MAX <= fabs(n))  // n ÷ d would overflow
            return r;
        return n / d;
    }

    void interpolateNewtonQuadratic(int stepCount)
    {
        // from Boost

        // trying to get a quadratic polynomial to fit the bracket
        // we also use the earlier interpolate ‘d’
        double A, B;
        A = safeDiv(fd - fb, d - b, DBL_MAX);
        B = safeDiv(fb - fa, b - a, DBL_MAX);
        A = safeDiv(A - B, d - a, 0);
        if (A == 0)  // failed
            return interpolateSecant();

        if (dontBracketARoot(A, fa))
            c = a;
        else
            c = b;

        // per SciPy we write the Newton method so that if one step falls in the bracket
        // and the next step falls out of it, the result of the earlier step is retained
        int step;
        for (step = 0; step < stepCount; ++step)
        {
            double c1 = c - safeDiv(fa + (B + A * (c - b)) * (c - a), B + A * (2 * c - a - b), 1 + c - a);
            if (outOfBracket(c1))
                break;
            c = c1;
        }
        if (step == 0)  // not even one step was successful
            interpolateSecant();
    }

    void interpolateCubic()
    {
        // from Boost

        double q11 = (d - e) * fd / (fe - fd);
        double q21 = (b - d) * fb / (fd - fb);
        double q31 = (a - b) * fa / (fb - fa);
        double d21 = (b - d) * fd / (fd - fb);
        double d31 = (a - b) * fb / (fb - fa);

        double q22 = (d21 - q11) * fb / (fe - fb);
        double q32 = (d31 - q21) * fa / (fd - fa);
        double d32 = (d31 - q21) * fd / (fd - fa);
        double q33 = (d32 - q22) * fa / (fe - fa);

        c = q31 + q32 + q33 + a;
    }

    bool cantDoCubic()
    {
        // from Boost; SciPy instead uses the AreClose params provided for the root

        static const double tol = eps * 32;
        return e != e ||  // fast isnan, true for first interpolation of second iteration
               fabs(fa - fb) < tol ||
               fabs(fa - fd) < tol ||
               fabs(fa - fe) < tol ||
               fabs(fb - fd) < tol ||
               fabs(fb - fe) < tol ||
               fabs(fd - fe) < tol;
    }

    bool evalRebracketAndCheckEnd()
    {
        // following tests from Boost
#if 0
        // the first test should not be required since we don't accept a tolerance below 4 × eps
        if (b - a < 4 * eps)
            interpolateBisection();
        else
#endif
        {
            double closestAllowedCNearA = a + 2 * eps * fabs(a);
            if (c < closestAllowedCNearA)
                c = closestAllowedCNearA;
            else
            {
                double closestAllowedCNearB = b - 2 * eps * fabs(b);
                if (c > closestAllowedCNearB)
                    c = closestAllowedCNearB;
            }
        }

        if (callFunctionAndCheckEnd(c, fc))
            return true;  // reached root or error

        // rebracket
        e = d; fe = fd;
        if (dontBracketARoot(fa, fc))
        {
            d = a; fd = fa;
            a = c; fa = fc;
        }
        else
        {
            d = b; fd = fb;
            b = c; fb = fc;
        }

        return checkBracketForEnd();
    }
};

extern "C"
double apsfindCustom(
    ApsFindInputFunction fn,
    void * otherInput,
    double intervalStart,
    double intervalEnd,
    ApsFindResultStatus * resultStatus,
    double absoluteTolerance,
    double relativeTolerance,
    int maximumIterations,
    int interpolationsPerIteration)
{
    int errorCode = APSFIND_NO_ERROR;
    if (!isfinite(intervalStart))
        errorCode |= APSFIND_INVALID_INTERVAL_START;
    if (!isfinite(intervalEnd))
        errorCode |= APSFIND_INVALID_INTERVAL_END;
    if (intervalStart >= intervalEnd)
        errorCode |= APSFIND_INVALID_INTERVAL;
    if (AreClose::invalidTolerance(absoluteTolerance))
        errorCode |= APSFIND_INVALID_ABSOLUTE_TOLERANCE;
    if (AreClose::invalidTolerance(relativeTolerance))
        errorCode |= APSFIND_INVALID_RELATIVE_TOLERANCE;
    if (maximumIterations < 1)
        errorCode |= APSFIND_INVALID_MAXIMUM_ITERATIONS;
    if (interpolationsPerIteration < 1)
        errorCode |= APSFIND_INVALID_INTERPOLATIONS_PER_ITERATION;

    if (errorCode != APSFIND_NO_ERROR)
    {
        if (resultStatus)
            resultStatus->errorCode = errorCode;
        return NAN;
    }

    return ApsFindSolver(fn, otherInput, intervalStart, intervalEnd, resultStatus, absoluteTolerance, relativeTolerance, maximumIterations, interpolationsPerIteration).result();
}

extern "C"
double apsfind(
    ApsFindInputFunction fn,
    void * otherInput,
    double intervalStart,
    double intervalEnd,
    ApsFindResultStatus * resultStatus)
{
    return apsfindCustom(fn, otherInput, intervalStart, intervalEnd, resultStatus, defaultAbsoluteTolerance, defaultRelativeTolerance, defaultMaximumIterations, defaultInterpolationsPerIteration);
}

typedef struct
{
    ApsFindUniFunction f;
    double target;
} _ApsFindUniParams;

static double _apsfindUniHelper(double guessInput, void * functionAndTarget)
{
    _ApsFindUniParams * p = static_cast<_ApsFindUniParams *>(functionAndTarget);
    return p->f(guessInput) - p->target;
}

extern "C"
double apsfindu(
    ApsFindUniFunction fn,
    double target,
    double intervalStart,
    double intervalEnd)
{
    _ApsFindUniParams params;
    params.f = fn;
    params.target = target;
    return apsfind(_apsfindUniHelper, (void *)(&params), intervalStart, intervalEnd, NULL);
}

extern "C"
void apsfindResultStatusPrint(FILE * file, ApsFindResultStatus stat, int precision)
{
    fprintf(file, "Iterations: %d, Function Calls: %d\n", stat.iterations, stat.functionCalls);
    fprintf(file, "Bracket: (%.*g, %.*g)\n", precision, stat.bracketStart, precision, stat.bracketEnd);
    if (stat.errorCode == 0)
        fprintf(file, "Error: none\n");
    // input errors, using & as they may be OR-ed together
    if (stat.errorCode & APSFIND_INVALID_INTERVAL_START)
        fprintf(file, "Error: Non-finite start of interval\n");
    if (stat.errorCode & APSFIND_INVALID_INTERVAL_END)
        fprintf(file, "Error: Non-finite end of interval\n");
    if (stat.errorCode & APSFIND_INVALID_INTERVAL)
        fprintf(file, "Error: Interval start should be less than interval end\n");
    if (stat.errorCode & APSFIND_INVALID_ABSOLUTE_TOLERANCE)
        fprintf(file, "Error: Invalid absolute tolerance, should be zero or finite and at least 4 × machine epsilon\n");
    if (stat.errorCode & APSFIND_INVALID_RELATIVE_TOLERANCE)
        fprintf(file, "Error: Invalid relative tolerance, should be zero or finite and at least 4 × machine epsilon\n");
    if (stat.errorCode & APSFIND_INVALID_MAXIMUM_ITERATIONS)
        fprintf(file, "Error: Maximum iterations should be at least 1\n");
    if (stat.errorCode & APSFIND_INVALID_INTERPOLATIONS_PER_ITERATION)
        fprintf(file, "Error: Interpolations per iteration should be at least 1\n");
    // execution errors, using == as they are mutually exclusive
    if (stat.errorCode == APSFIND_INVALID_FUNCTION_VALUE)
        fprintf(file, "Error: Non-finite function value encountered\n");
    if (stat.errorCode == APSFIND_INTERVAL_DOES_NOT_BRACKET_A_ROOT)
        fprintf(file, "Error: Interval does not bracket a root\n");
    if (stat.errorCode == APSFIND_MAXIMUM_ITERATIONS_REACHED)
        fprintf(file, "Error: Maximum iterations reached\n");
}
