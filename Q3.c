#include <complex.h>
#include <stdio.h>
#include <time.h>

/*
  Note for self: Understanding spectral leakage: When we perform an FFT, we are multiplying each signal window
  with a corresponding sin and cosine wave of n frequency. However, our frequency resolution is limited
  by our sampling rate/sample length. If this resolution is not fine enough to correspond with exactly with
  the nth sin/cosine signal present it will not be exactly correlated, and thus, it will also correlated a bit
  with neighboring frequency samples leading to spectral leakage. Furthermore, due to nature of real signals, even if there is correlation
  between sin/cosine wave of nth frequency, it is possible that in the real signal, there is only a
  portion of that sin/cosine way (DFT assumes it will complete/periodic). This gap between signals when making correlations leads to introduction of noisy
  high frequency components (how exactly? not so sure about this) and consequently cause spectral leakage.

  Two ways to address this:
  1. Make sure our signal sampling resolution is fine enough to capture closely/accurately present frequencies
  2. Use windowing techniques to avoid gaps when performing FFT/DFT.
*/

// Global definition for pi
const double pi = 3.141592654;

// Function declarations
double normalizeRad(double x);
double myCos(double x);
double *convolution(double *x, int sampleSizeX, double *kernel, int kernelSize);
double *generateHanningWindow(int N);

int main() {
    double x[] = {1, 1, 1, 0, 0, 0 ,1};
    int sampleSizeX = sizeof(x)/sizeof(x[0]);
    //double kernel[] = {0.5, 0.2, 0.3, 0.5, 0.2, 0.3, 0.5, 0.2, 0.3 };
    //int kernelSize = sizeof(kernel)/sizeof(kernel[0]);
    int kernelSize = 3;
    double *kernel = generateHanningWindow(kernelSize);

    double *result = convolution(x, sampleSizeX, kernel, kernelSize);

    for (int i = 0; i < sampleSizeX + kernelSize - 1; i++) {
      printf("%f\n", result[i]);
    }

    free(result);
    //free(kernel);
}

double *generateHanningWindow(int N) {
  /*
    This function generates a 1D hanning window of size N
  */
  double *result = malloc((N)*sizeof(double));
  for (int i = 0; i < N; i++) {
    result[i] = 0.5*(1-myCos(2*pi*i/N));
  }

  return result;
}

double *convolution(double *x, int sampleSizeX, double *kernel, int kernelSize) {
  /* 
    This function performs convolution between a signal/array and a passed kernel. 
  */

  double *result = malloc((sampleSizeX + kernelSize - 1)*sizeof(double)); // Array to store results

  /* 
    Iterate over each sample, summing the possible products of the signal and
    reversed kernel at that index as it slides over the signal .
  */
  for (int i = 0; i < sampleSizeX + kernelSize - 1; i++) {
      result[i] = 0;
      for (int j = 0; j < kernelSize; j++) {
          if (i-j >= 0 && i-j < sampleSizeX) { 
              /* 
                The condition ensures that kernel index combinations that do not correspond
                to valid indexes on the signal are ignored.
              */
              //printf("%f, %f, %i, %i, %f\n", x[i-j], kernel[j], i-j, j, x[i-j] * kernel[j]);
              result[i] += x[i-j] * kernel[j];
          }
      }
      //printf("\n", result[i]);
  }

  return result;
}


double myCos(double x) {
  /*
    This function uses a Taylor Series expansion
    to approximate cos(x)
  */

  // Normalizing value between [-pi, pi] since
  // that is the range in which our approximation lies
  x = normalizeRad(x);
  int approxOrder = 7;  // Approximation order
  int positiveTerm = 0; // Acts a boolean to alternate between +ve & -ve terms
                        // in the approximation

  // Keeping track of these terms to avoid redundant computation
  int factorialTerm = 1;
  double powerTerm = 1;

  double result = 1;

  // Keeping track of these terms to avoid redundant computation
  for (int i = 2; i < approxOrder * 2; i += 2) {
    if (positiveTerm) {
      powerTerm = x * x * powerTerm;
      factorialTerm = i * (i - 1) * factorialTerm;
      result += (powerTerm / factorialTerm);
      positiveTerm = 0;
    } else {
      powerTerm = x * x * powerTerm;
      factorialTerm = i * (i - 1) * factorialTerm;
      result -= (powerTerm / factorialTerm);
      positiveTerm = 1;
    }
  }

  return result;
}

double normalizeRad(double x) {
  /*
    Taking input a double and normalizes it to range [-pi, pi]
  */

  if (x > pi) {
    while (x > pi) {
      x -= 2 * pi; // If we have positive double greater than pi, keep
                   // subtracting until it yes than pi
    }
  } else if (x < -pi) {
    while (x < -pi) {
      x += 2 * pi;
    }
  }

  return x;
}