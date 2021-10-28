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

#ifdef __cplusplus
extern "C"
{
#endif

// type declarations

enum BisectionErrorCode
{
    // errors in function input, may be OR-ed together
    BISECTION_NO_ERROR                             =   0,
    BISECTION_INVALID_INTERVAL_START               =   1,
    BISECTION_INVALID_INTERVAL_END                 =   2,
    BISECTION_INVALID_INTERVAL                     =   4,
    BISECTION_INVALID_ABSOLUTE_TOLERANCE           =   8,
    BISECTION_INVALID_RELATIVE_TOLERANCE           =  16,
    BISECTION_INVALID_MAXIMUM_ITERATIONS           =  32,
    // errors in function execution, mutually exclusive
    BISECTION_INVALID_FUNCTION_VALUE               = 128,
    BISECTION_INTERVAL_DOES_NOT_BRACKET_A_ROOT     = 129,
    BISECTION_MAXIMUM_ITERATIONS_REACHED           = 130
};

typedef struct
{
    int iterations, functionCalls;
    double bracketStart, bracketEnd;
    int errorCode;
} BisectionResultStatus;

typedef double(*BisectionInputFunction)(double, void *);

// function declarations

double bisectionCustom(
    BisectionInputFunction function,
    void * otherInput,
    double intervalStart,
    double intervalEnd,
    BisectionResultStatus * resultStatus,
    double absoluteTolerance,
    double relativeTolerance,
    int maximumIterations);

double bisection(
    BisectionInputFunction function,
    void * otherInput,
    double intervalStart,
    double intervalEnd,
    BisectionResultStatus * resultStatus);

#ifdef __cplusplus
} // extern "C"
#endif
