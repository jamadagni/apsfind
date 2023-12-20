// jaya jaya shankara hara hara shankara

#include "../apsfind.h"
#include <math.h>
#include <stdio.h>

/*
The formula of the superellipse is:
(x / a) ** n + (y / b) ** n = 1

The area of the superellipse is:

4 * a * b * pow(tgamma(1 + 1 / n), 2) / tgamma(1 + 2 / n)

For a supercircle a = b and for the unit supercircle a = b = 1.

The area of a unit supercircle is hence as below.
*/

double areaOfUnitSuperCircle(double n)
{
    return 4 * pow(tgamma(1 + 1 / n), 2) / tgamma(1 + 2 / n);
}

int main()
{
    double circleVal, meanVal;

    // Consider the circle of unit radius (and centered at the origin).
    // Consider its bounding square.

    // Area of the circle is 3.1416~ (pi)
    // As the circle is a supercircle of degree 2,
    // let's obtain the area that way just to illustrate:
    circleVal = areaOfUnitSuperCircle(2);
    printf("Area of the unit circle     = %f (pi)\n", circleVal);

    // Now area of the above square is 4.
    printf("Area of its bounding square = 4\n");

    // Let's find out the mean between the area of the circle and the square
    meanVal = (circleVal + 4) / 2;
    printf("The mean of the above       = %f\n", meanVal);

    printf("The degree of the unit supercircle with above as area:\n");
    printf("%f\n", apsfindu(areaOfUnitSuperCircle, meanVal, 2, 10));

    printf("Note that this is close to but not equal to pi!\n");

    // REF:
    // Gardner, Martin, 1977, “Piet Hein's Superellipse”, Mathematical Carnival, Vintage Press, New York.
    // The idea of a superellipse with mean area between a normal ellipse and its bounding rectangle is
    // mentioned on p 250, as is the observation that the degree is close to pi. However the article says
    // that the required degree “is a trifle greater than 3.17” whereas we find it as above to be ~3.162.
}
