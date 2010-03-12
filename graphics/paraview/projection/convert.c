#include <math.h>

void CartesianToSpherical(float x, float y, float z, float* rho, float* phi, float* theta)
{
        *rho = sqrtf( (x*x) + (y*y) + (z*z));
        *theta = atan2f(y, x);
        *phi = acosf(z/(*rho));
}

void SphericalToCartesian(float rho, float phi, float theta, float* x, float* y, float* z)
{
        *x = rho* sinf(phi) * cosf(theta);
        *y = rho* sinf(phi) * sinf(theta);
        *z = rho* cosf(phi);
}

