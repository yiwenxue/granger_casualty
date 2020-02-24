#ifndef CASUALTY_NOISE
#define CASUALTY_NOISE
#include <ctime>
#include <fftw3.h>
#include <mathematics.h>

double * pinkNoise(double exponent, int size);

int GenpinkNoise(double *noise, double exponent, int size);

int test_noise(void);



#endif
