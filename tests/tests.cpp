// Copyright © 2006, John Maddock
// Copyright © 2021, Shriramana Sharma

// Test function implements the same test cases as used by
// "Algorithm 748: Enclosing Zeros of Continuous Functions"
// by G.E. Alefeld, F.A. Potra and Yixun Shi.

// Following implementation adapted from:
// https://github.com/boostorg/math/blob/master/test/test_toms748_solve.cpp

// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0

#ifdef TOMS748
#include "toms748.h"
#endif

#ifdef BISECTION
#include "bisection.h"
#endif

#include <cmath>
#include <cstdio>

using namespace std;

struct TestParams {
    TestParams(int testID): t{testID} {}
    TestParams(int testID, int ival): t{testID}, i{ival} {}
    TestParams(int testID, double jval, double kval): t{testID}, j{jval}, k{kval} {}
    int t, i;
    double j, k;
};

double test(double x, void * otherInput)
{
    auto tp = static_cast<TestParams *>(otherInput);
    int t = tp->t, i = tp->i;
    double j = tp->j, k = tp->k;

    switch(t)
    {
    case 1:
        return sin(x) - x / 2;
    case 2:
        {
            double r = 0;
            for(int i = 1; i <= 20; ++i)
            {
                double p = (2 * i - 5);
                double q = x - i * i;
                r += p * p / (q * q * q);
            }
            r *= -2;
            return r;
        }
    case 3:
        return j * x * exp(k * x);
    case 4:
        return pow(x, k) - j;
    case 5:
        return sin(x) - 0.5;
    case 6:
        return 2 * x * exp(-double(i)) - 2 * exp(-i * x) + 1;
    case 7:
        return (1 + (1 - i) * (1 - i)) * x - (1 - i * x) * (1 - i * x);
    case 8:
        return x * x - pow(1 - x, i);
    case 9:
        return (1 + (1 - i) * (1 - i) * (1 - i) * (1 - i)) * x - (1 - i * x) * (1 - i * x) * (1 - i * x) * (1 - i * x);
    case 10:
        return exp(-i * x) * (x - 1) + pow(x, double(i));
    case 11:
        return (i * x - 1) / ((i - 1) * x);
    case 12:
        return pow(x, double(1)/i) - pow(double(i), double(1)/i);
    case 13:
        return x == 0 ? 0 : x / exp(1 / (x * x));
    case 14:
        return x >= 0 ? (double(i) / 20) * (x / 1.5f + sin(x) - 1) : -double(i) / 20;
    case 15:
        {
            double d = 2e-3 / (1 + i);
            if (x > d)
                return exp(1.0) - 1.859;
            else if(x > 0)
                return exp((i + 1) * x * 1000 / 2) - 1.859;
            return -0.859;
        }
    }
    return 0;
}

// Following test result values for comparison adapted from SciPy:
// https://github.com/scipy/scipy/blob/master/scipy/optimize/_tstutils.py

double test_results_from_scipy[] = {
    1.895494267033981, 3.0229153472730568, 6.6837535608080785,
    11.238701655002211, 19.67600008062341, 29.828227326504756,
    41.906116195289414, 55.95359580014309, 71.9856655865878,
    90.00886853916667, 110.0265327483302, 0,
    0, 0, 0.668740304976422,
    0.76472449133173, 0.8177654339579425, 0.8513399225207846,
    0.8744852722211679, 1, 1,
    1, 1, 1,
    1, 1, 1,
    1, 0.5235987755982988, 0.4224777096412367,
    0.3066994104832037, 0.22370545765466296, 0.17171914751950837,
    0.13825715505682407, 0.03465735902085385, 0.01732867951399863,
    0.011552453009332421, 0.008664339756999316, 0.006931471805599454,
    0.0384025518406219, 0.0099000099980005, 0.0024937500390620117,
    0.5, 0.34595481584824206, 0.24512233375330722,
    0.19554762353656563, 0.16492095727644096, 0.2755080409994844,
    0.1377540204997422, 0.010305283778156442, 0.0036171081789040634,
    0.0004108729184963954, 2.598957589290763e-05, 7.668595122185337e-06,
    0.401058137541547, 0.5161535187579336, 0.5395222269084158,
    0.5481822943406552, 0.5527046666784878, 0.5,
    0.2, 0.06666666666666667, 0.05,
    2, 3, 4,
    5, 6, 7,
    9, 11, 13,
    15, 17, 19,
    21, 23, 25,
    27, 29, 31,
    33, 0, 0.6238065189616124,
    0.6238065189616124, 0.6238065189616124, 0.6238065189616124,
    0.6238065189616124, 0.6238065189616124, 0.6238065189616124,
    0.6238065189616124, 0.6238065189616124, 0.6238065189616124,
    0.6238065189616124, 0.6238065189616124, 0.6238065189616124,
    0.6238065189616124, 0.6238065189616124, 0.6238065189616124,
    0.6238065189616124, 0.6238065189616124, 0.6238065189616124,
    0.6238065189616124, 0.6238065189616124, 0.6238065189616124,
    0.6238065189616124, 0.6238065189616124, 0.6238065189616124,
    0.6238065189616124, 0.6238065189616124, 0.6238065189616124,
    0.6238065189616124, 0.6238065189616124, 0.6238065189616124,
    0.6238065189616124, 0.6238065189616124, 0.6238065189616124,
    0.6238065189616124, 0.6238065189616124, 0.6238065189616124,
    0.6238065189616124, 0.6238065189616124, 0.6238065189616124,
    5.905130559421972e-05, 5.6367155339937e-05, 5.391640945559192e-05,
    5.166989239494225e-05, 4.960309669914456e-05, 4.7695285287638995e-05,
    4.5928793239948666e-05, 4.4288479195664784e-05, 4.276129025788324e-05,
    4.13359139159538e-05, 4.000249733801981e-05, 3.875241929620669e-05,
    3.757810355995799e-05, 3.6472865219959236e-05, 3.543078335653183e-05,
    3.44465949299615e-05, 3.351560587780037e-05, 3.263361624943721e-05,
    3.1796856858426e-05, 3.1001935436965345e-05, 3.0245790670210097e-05,
    1.2277994232461523e-05, 6.169539390440866e-06, 4.119858529829282e-06,
    3.092462387727217e-06, 2.475204426105018e-06, 2.063356767851271e-06,
    1.7690120078154265e-06, 1.5481615698859102e-06, 1.3763345366022351e-06,
    1.238838578899714e-06};
int test_results_from_scipy_index = 0;


void run_test_helper(TestParams &tp, double a, double b)
{
#ifdef TOMS748
    Toms748ResultStatus result_status;
    double r = toms748(test, static_cast<void *>(&tp), a, b, &result_status);
#endif

#ifdef BISECTION
    BisectionResultStatus result_status;
    double r = bisection(test, static_cast<void *>(&tp), a, b, &result_status);
#endif

    double expectedResult = test_results_from_scipy[test_results_from_scipy_index];
    printf("Function %2d ; ", tp.t);
    printf("result = %10.6e %10.6f ; ", r, r);
    printf("deviation = %6.2f%% ; ", expectedResult == 0 ? fabs(r - expectedResult) * 100 : fabs((r - expectedResult) / expectedResult) * 100);
    printf("iters = %2d ; ", result_status.iterations);
    printf("funcalls = %2d\n", result_status.functionCalls);
    test_results_from_scipy_index += 1;
}

void run_test(int id, double a, double b)
{
    auto tp = TestParams(id);
    run_test_helper(tp, a, b);
}

void run_test(int id, double a, double b, int p)
{
    auto tp = TestParams(id, p);
    run_test_helper(tp, a, b);
}

void run_test(int id, double a, double b, double p1, double p2)
{
    auto tp = TestParams(id, p1, p2);
    run_test_helper(tp, a, b);
}

int main()
{
    run_test(1, 3.14/2, 3.14);

    for(int i = 1; i <= 10; i += 1)
    {
        run_test(2, i*i + 1e-9, (i+1)*(i+1) - 1e-9);
    }

    run_test(3, -9.0, 31.0, -40.0, -1.0);
    run_test(3, -9.0, 31.0, -100.0, -2.0);
    run_test(3, -9.0, 31.0, -200.0, -3.0);

    for(int n = 4; n <= 12; n += 2)
    {
        run_test(4, 0.0, 5.0, 0.2, double(n));
    }
    for(int n = 4; n <= 12; n += 2)
    {
        run_test(4, 0.0, 5.0, 1.0, double(n));
    }
    for(int n = 8; n <= 14; n += 2)
    {
        run_test(4, -0.95, 4.05, 1.0, double(n));
    }
    run_test(5, 0.0, 1.5);
    for(int n = 1; n <= 5; ++n)
    {
        run_test(6, 0.0, 1.0, n);
    }
    for(int n = 20; n <= 100; n += 20)
    {
        run_test(6, 0.0, 1.0, n);
    }
    run_test(7, 0.0, 1.0, 5);
    run_test(7, 0.0, 1.0, 10);
    run_test(7, 0.0, 1.0, 20);
    run_test(8, 0.0, 1.0, 2);
    run_test(8, 0.0, 1.0, 5);
    run_test(8, 0.0, 1.0, 10);
    run_test(8, 0.0, 1.0, 15);
    run_test(8, 0.0, 1.0, 20);
    run_test(9, 0.0, 1.0, 1);
    run_test(9, 0.0, 1.0, 2);
    run_test(9, 0.0, 1.0, 4);
    run_test(9, 0.0, 1.0, 5);
    run_test(9, 0.0, 1.0, 8);
    run_test(9, 0.0, 1.0, 15);
    run_test(9, 0.0, 1.0, 20);
    run_test(10, 0.0, 1.0, 1);
    run_test(10, 0.0, 1.0, 5);
    run_test(10, 0.0, 1.0, 10);
    run_test(10, 0.0, 1.0, 15);
    run_test(10, 0.0, 1.0, 20);

    run_test(11, 0.01, 1.0, 2);
    run_test(11, 0.01, 1.0, 5);
    run_test(11, 0.01, 1.0, 15);
    run_test(11, 0.01, 1.0, 20);

    for(int n = 2; n <= 6; ++n)
        run_test(12, 1.0, 100.0, n);
    for(int n = 7; n <= 33; n+=2)
        run_test(12, 1.0, 100.0, n);

    run_test(13, -1.0, 4.0);

    for(int n = 1; n <= 40; ++n)
        run_test(14, -1e4, 3.14/2, n);

    for(int n = 20; n <= 40; ++n)
        run_test(15, -1e4, 1e-4, n);

    for(int n = 100; n <= 1000; n+=100)
        run_test(15, -1e4, 1e-4, n);
}
