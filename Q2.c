#include <complex.h>
#include <stdio.h>
#include <time.h>

// Global definition for pi
const double pi = 3.141592654;

// Function declarations
double normalizeRad(double x);
double mySin(double x);
double myCos(double x);
double complex euler(double complex x);
double complex *DFT(double complex x[], int sampleSize);
double complex *FFT(double complex x[], int sampleSize, int binSize);
double complex *FFT_Recursion(double complex x[], int sampleSize);

int main() {

  // Using clock to keep track of how long it takes to compute normal DFT vs FFT
  clock_t start, end;

  // double complex x[] = {2, 1, -1, 5, 0, 3, 0, -4, 2, 1, -1, 5, 0, 3, 0, -4,
  // 2, 1, -1, 5, 0, 3, 0, -4, 2, 1, -1, 5, 0, 3, 0, -4};
  double complex x[] = {2, 1, -1, 5, 0, 3, 0, -4};
  int sampleSize =
      sizeof(x) / sizeof(x[0]); // Calculating sample or input signal size
  start = clock();
  double complex *resultArray = DFT(x, sampleSize); // Performing DFT
  end = clock();

  // printf("CPU time used: %f seconds\n",
  //        ((double)(end - start)) /
  //            CLOCKS_PER_SEC); // Printing time taken to compute
  // printf("\n");

  // double complex y[] = {2, 1, -1, 5, 0, 3, 0, -4, 2, 1, -1, 5, 0, 3, 0, -4,
  // 2, 1, -1, 5, 0, 3, 0, -4, 2, 1, -1, 5, 0, 3, 0, -4};
  double complex y[] = {2, 1, -1, 5, 0, 3, 0, -4};
  sampleSize = sizeof(y) / sizeof(y[0]);
  int binSize = 8;
  start = clock();
  FFT(y, sampleSize, 8);
  end = clock();

  // printf("CPU time used: %f seconds\n",
  //        ((double)(end - start)) / CLOCKS_PER_SEC);
}

double complex *FFT(double complex x[], int sampleSize, int binSize) {
  /*
    This function calculates the FFT of a one-dimenional input signal. It takes
    input the sampleSize, which is the size of the signal and the binSize, which
    reflects the sample rate at which we want to perform FFT on the signal. It
    is worth noting that this function works on the assumpton that the signal is
    of length 2^n (a way to fix this is to add zero-padding to the signal to
    force it into an order 2^n).

    Notes for self: If the rate at which the signal has been
    sampled and the rate at which the FFT has been calculated are not the same,
    this can lead to spectral leakage, which corre
  */

  // Calculating bins needed for each signal sub-sample
  int bins = sampleSize / binSize;
  double complex *result[bins];

  // Calculating FFT for each sub-sample
  for (int i = 0; i < bins; i++) {
    result[i] = FFT_Recursion(x + i * binSize, binSize);
    printf("[");
    for (int j = 0; j < binSize; j++) {
      printf("%.2f + %.2fi, ", creal(result[i][j]), cimag(result[i][j]));
    }
    printf("]\n");
  }

  return result;
}

double complex *FFT_Recursion(double complex x[], int sampleSize) {
  /*
    This function represents a recursive function that allows us to calculate
    the DFT using a time-complexity of n*log(n) instead of n^2, which
    constitutes a simplifed approach.

    Notes for self: the FFT achieves such efficiency by exploiting symmetry in
    Fourier Matrix (DFT), which leads to common twiddle factors that allow us to
    relate the odd and even components while calculating the FFT. This allowing
    us to break the each DFT step into a smaller requiring step (N^2/2)
    computations. Recursively, this leads to n log n computations.

    The magic here happens because of symmetries in the odd and even components, where for
    each sub-calcuation we see that beyond e^-j*pi*k*n/N/2 for k > N/2, the result is symmetric.
    This allows use to break the step in N^2/4 (for even) & N^2/4 (for odd) computations at each step. 
  */
  if (sampleSize <= 1) {
    return x; // Base case, can no longer recursively break
  } else {

    // Dynamically allocated memory for odd and even samples so that our array's
    // aren't lost once the scope of a function call ends.
    double complex *even = malloc(sampleSize / 2 * sizeof(double complex));
    double complex *odd = malloc(sampleSize / 2 * sizeof(double complex));

    // Splitting our sample into separate odd and even components
    for (int i = 0; i < sampleSize / 2; i++) {
      even[i] = x[i * 2];
      odd[i] = x[(i * 2) + 1];
    }

    // Applying fft recursively to further divide the problem
    FFT_Recursion(even, sampleSize / 2);
    FFT_Recursion(odd, sampleSize / 2);

    // Once we hit the base case, we go back upwards from here
    for (int i = 0; i < sampleSize / 2; i++) {
      // Calculating the twiddle factor that allows us to relate the odd and
      // even components
      double complex twiddlePlusOdd = euler((2 * pi * i) / sampleSize) * odd[i];
      x[i] =
          even[i] + twiddlePlusOdd; // Re-joining even and odd component for first half
      x[i + sampleSize / 2] =
          even[i] - twiddlePlusOdd; // Re-joining even and odd component for second half
    }

    free(even); // Freeing the dynamically allocated memory
    free(odd);
  }
}

double complex *DFT(double complex x[], int sampleSize) {
  /*
    This function implements the raw DFT which is not optimized.

    Notes for self: the DFT is successful at extracting frequencies because it
    multiplies the signal with cosine and sin terms are different frequencies capturing
    the correlation as an amplitude. Thus, the frequency components that match have a higher
    amplitude and contribute largely to the final result. 

    Q: How do we ensure that we capture or compare/correlate with all possible frequencies?
    There is possiblility that if our sample size is small, our FFT fails to capture frequency
    components present in the signal?
    Nyquist Theorm: sample at twice the highest frequency present in the signal. The highest
    frequency present in the signal will correspond directly to the frequency at which the signal
    was sampled.  
  */
  for (int i = 0; i < sampleSize; i++) {
    double complex Xk = 0 + 0 * I;
    for (int j = 0; j < sampleSize; j++) {
      Xk += x[j] * euler((2 * pi * i * j) / sampleSize);
    }
    x[i] = Xk;
  }

  return x;
}

double complex euler(double complex x) {
  /*
    Uses the euler identity for calculating the value of complex exponentials
  */
  return myCos(creall(x)) + (I * mySin(creall(x)));
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

double mySin(double x) {
  /*
    This function uses a Taylor Series expansion
    to approximate cos(x)
  */

  // Normalizing value between [-pi, pi] since
  // that is the range in which our approximation lies
  x = normalizeRad(x);
  int approxOrder = 6;  // Approximation order
  int positiveTerm = 0; // Acts a boolean to alternate between +ve & -ve terms
                        // in the approximation

  // Keeping track of these terms to avoid redundant computation
  int factorialTerm = 1;
  double powerTerm = x;
  double result = x;

  for (int i = 3; i < approxOrder * 2; i += 2) {
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