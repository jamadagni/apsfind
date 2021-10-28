// bisection
// =========
//
// Bisection algorithm to find a root of a function in a given interval.
//
// This is a C++ implementation with a C-compatible interface.
// Copyright © 2021, Shriramana Sharma, samjnaa-at-gmail-dot-com
//
// Use, modification and distribution are permitted subject to the
// "BSD-2-Clause"-type license stated in the accompanying file LICENSE.txt

#include "bisection.h"
#include "areclose.hpp"
#include <cmath>
#include <cfloat>
using namespace std;


// default values below from SciPy
static const int
    defaultMaximumIterations = 100;
static const double
    NaN = nan(""),
    defaultRelativeTolerance = DBL_EPSILON * 4,
    defaultAbsoluteTolerance = 2e-12;


// The idea of the generic pointer for other fixed input is from:
// SciPy Cython: https://docs.scipy.org/doc/scipy/reference/optimize.cython_optimize.html
//
// The ResultStatus struct combines SciPy's RootResults class and Boost's returning the final bracket

struct BisectionSolver
{
private:
    BisectionInputFunction f;
    void * args;
    double a, fa, b, fb, c, fc;
    BisectionResultStatus * stat;
    AreClose areClose;
    int maxiter, remiter, fncount;
    double root = NaN;

public:
    BisectionSolver(
        BisectionInputFunction function,
        void * otherInput,
        double intervalStart,
        double intervalEnd,
        BisectionResultStatus * resultStatus,
        double absoluteTolerance,
        double relativeTolerance,
        int maximumIterations):
        f{function},
        args{otherInput},
        a{intervalStart},
        b{intervalEnd},
        stat{resultStatus},
        areClose{absoluteTolerance, relativeTolerance},
        maxiter{maximumIterations},
        remiter{maximumIterations},
        fncount{0}
    {
        if (checkBracketForEnd())  // maybe already input is so
            return;
        if (callFunctionAndCheckEnd(a, fa))
            return;
        if (callFunctionAndCheckEnd(b, fb))
            return;
        if (dontBracketARoot(fa, fb))
        {
            finalize(NaN, BISECTION_INTERVAL_DOES_NOT_BRACKET_A_ROOT, a, b);
            return;
        }

        // iteration loop

        for (; remiter; --remiter)
        {
            interpolateBisection();
            if (evalRebracketAndCheckEnd())
                return;
        }

        finalize(NaN, BISECTION_MAXIMUM_ITERATIONS_REACHED, a, b);
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
        finalize(c, BISECTION_NO_ERROR, a, b);
        return true;
    }

    bool callFunctionAndCheckEnd(double x, double & fx)
    {
        fx = f(x, args);
        fncount += 1;
        if (fx == 0)
        {
            finalize(x, BISECTION_NO_ERROR, x, x);
            return true;
        }
        if (!isfinite(fx))
        {
            finalize(NaN, BISECTION_INVALID_FUNCTION_VALUE, a, b);
            return true;
        }
        return false;
    }

    inline static bool dontBracketARoot(double p, double q)
    {
        return p * q > 0;
    }

    inline bool outOfBracket(double val)
    {
        return val <= a || val >= b;
    }

    inline void interpolateBisection()
    {
        // from Boost; SciPy just uses (a + b) / 2
        c = a + (b - a) / 2;
    }

    bool evalRebracketAndCheckEnd()
    {
        if (callFunctionAndCheckEnd(c, fc))
            return true;  // reached root or error

        // rebracket
        if (dontBracketARoot(fa, fc))
        {
            a = c; fa = fc;
        }
        else
        {
            b = c; fb = fc;
        }

        return checkBracketForEnd();
    }
};

extern "C"
double bisectionCustom(
    BisectionInputFunction function,
    void * otherInput,
    double intervalStart,
    double intervalEnd,
    BisectionResultStatus * resultStatus,
    double absoluteTolerance,
    double relativeTolerance,
    int maximumIterations)
{
    int errorCode = BISECTION_NO_ERROR;
    if (!isfinite(intervalStart))
        errorCode |= BISECTION_INVALID_INTERVAL_START;
    if (!isfinite(intervalEnd))
        errorCode |= BISECTION_INVALID_INTERVAL_END;
    if (intervalStart >= intervalEnd)
        errorCode |= BISECTION_INVALID_INTERVAL;
    if (AreClose::invalidTolerance(absoluteTolerance))
        errorCode |= BISECTION_INVALID_ABSOLUTE_TOLERANCE;
    if (AreClose::invalidTolerance(relativeTolerance))
        errorCode |= BISECTION_INVALID_RELATIVE_TOLERANCE;
    if (maximumIterations < 1)
        errorCode |= BISECTION_INVALID_MAXIMUM_ITERATIONS;

    if (errorCode != BISECTION_NO_ERROR)
    {
        if (resultStatus)
            resultStatus->errorCode = errorCode;
        return NaN;
    }

    return BisectionSolver(function, otherInput, intervalStart, intervalEnd, resultStatus, absoluteTolerance, relativeTolerance, maximumIterations).result();
}

extern "C"
double bisection(
    BisectionInputFunction function,
    void * otherInput,
    double intervalStart,
    double intervalEnd,
    BisectionResultStatus * resultStatus)
{
    return bisectionCustom(function, otherInput, intervalStart, intervalEnd, resultStatus, defaultAbsoluteTolerance, defaultRelativeTolerance, defaultMaximumIterations);
}
