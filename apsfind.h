// apsfind
// =======
//
// The algorithm of Alefeld, Potra and Shi to find a root of a function in
// a given interval using inverse cubic and “Newton-quadratic” interpolations.
// Published as TOMS 748, see reference in source.
//
// This is a C++ implementation with a C-compatible interface and an optional C++ helper.
// Copyright © 2021, Shriramana Sharma, samjnaa-at-gmail-dot-com
//
// Use, modification and distribution are permitted subject to the
// "BSD-2-Clause"-type license stated in the accompanying file LICENSE.txt

#ifdef __cplusplus
extern "C"
{
#endif

// type declarations

enum ApsFindErrorCode
{
    // errors in function input, may be OR-ed together
    APSFIND_NO_ERROR                             =   0,
    APSFIND_INVALID_INTERVAL_START               =   1,
    APSFIND_INVALID_INTERVAL_END                 =   2,
    APSFIND_INVALID_INTERVAL                     =   4,
    APSFIND_INVALID_ABSOLUTE_TOLERANCE           =   8,
    APSFIND_INVALID_RELATIVE_TOLERANCE           =  16,
    APSFIND_INVALID_MAXIMUM_ITERATIONS           =  32,
    APSFIND_INVALID_INTERPOLATIONS_PER_ITERATION =  64,
    // errors in function execution, mutually exclusive
    APSFIND_INVALID_FUNCTION_VALUE               = 128,
    APSFIND_INTERVAL_DOES_NOT_BRACKET_A_ROOT     = 256,
    APSFIND_MAXIMUM_ITERATIONS_REACHED           = 512
};

typedef struct
{
    int iterations, functionCalls;
    double bracketStart, bracketEnd;
    int errorCode;
} ApsFindResultStatus;

typedef double(*ApsFindInputFunction)(double, void *);

typedef double(*ApsFindUniFunction)(double);

// function declarations

double apsfindCustom(
    ApsFindInputFunction function,
    void * otherInput,
    double intervalStart,
    double intervalEnd,
    ApsFindResultStatus * resultStatus,
    double absoluteTolerance,
    double relativeTolerance,
    int maximumIterations,
    int interpolationsPerIteration);

double apsfind(
    ApsFindInputFunction function,
    void * otherInput,
    double intervalStart,
    double intervalEnd,
    ApsFindResultStatus * resultStatus);

double apsfindu(
    ApsFindUniFunction function,
    double target,
    double intervalStart,
    double intervalEnd);

#include <stdio.h>
void apsfindResultStatusPrint(FILE * file, ApsFindResultStatus stat, int prec);

#ifdef __cplusplus
} // extern "C"

template<typename Functor>
double _apsfindWrapperToFunctor(double var, void * fixed)
{
    return static_cast<Functor *>(fixed)->operator()(var);
}

template<typename Functor>
double apsfind(
    Functor fr,
    double intervalStart,
    double intervalEnd,
    ApsFindResultStatus * stat = nullptr)
{
    return apsfind(&_apsfindWrapperToFunctor<Functor>, static_cast<void *>(&fr), intervalStart, intervalEnd, stat);
}
#endif
