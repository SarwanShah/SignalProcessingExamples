#include <complex.h>

// Function declarations
float normalizeRad(float x);
float mySin(float x);
float myCos(float x);
long double myExp(long double x);

// Global definition for pi
float pi = 3.141592654;


int main() {
  float x = (-5   * pi) / 2;
  float sinResult = mySin(x);
  float cosResult = myCos(x);
  float expResult = myExp(5);
  printf("Sin(x): %f\n", sinResult);
  printf("Cos(x): %f\n", cosResult);
  printf("Exp(x): %f\n", expResult);
}

long double myExp(long double x) {
  /* 
    This function uses a Taylor Series
    approximation to estimate e^x
  */

  int approxOrder = 15; // The order of the approximation (Higher gives better approximation).

  // Keeping track of the factorial and power term so we don't have to compute
  // it from scratch for following terms and can use previous results
  long int factorialTerm = 1; 
  long double powerTerm = x;

  // Computing first two terms
  long double result = 1 + x;

  // Using a loop to for higher order terms
  for (int i = 2; i < approxOrder; i++) {
    powerTerm = x*powerTerm;
    factorialTerm = i*factorialTerm;
    result += powerTerm / factorialTerm;
  }

  return result;
}

float myCos(float x) {
  /* 
    This function uses a Taylor Series expansion
    to approximate cos(x)
  */

  // Normalizing value between [-pi, pi] since
  // that is the range in which our approximation lies
  x = normalizeRad(x);
  int approxOrder = 7; // Approximation order
  int positiveTerm = 0; // Acts a boolean to alternate between +ve & -ve terms in the approximation

  // Keeping track of these terms to avoid redundant computation
  int factorialTerm = 1; 
  float powerTerm = 1;

  float result = 1;

  // Keeping track of these terms to avoid redundant computation
  for (int i = 2; i < approxOrder * 2; i += 2) {
    if (positiveTerm) {
      powerTerm = x*x*powerTerm; 
      factorialTerm = i*(i-1)*factorialTerm;   
      result += (powerTerm / factorialTerm);
      positiveTerm = 0;
    } else {
      powerTerm = x*x*powerTerm; 
      factorialTerm = i*(i-1)*factorialTerm;   
      result -= (powerTerm / factorialTerm);    
      positiveTerm = 1;
    }
  }

  return result;
}

float mySin(float x) {
  /* 
    This function uses a Taylor Series expansion
    to approximate cos(x)
  */

  // Normalizing value between [-pi, pi] since
  // that is the range in which our approximation lies
  x = normalizeRad(x);
  int approxOrder = 6; // Approximation order
  int positiveTerm = 0; // Acts a boolean to alternate between +ve & -ve terms in the approximation

  // Keeping track of these terms to avoid redundant computation
  int factorialTerm = 1;
  float powerTerm = x;
  float result = x;

  for (int i = 3; i < approxOrder * 2; i += 2) {
    if (positiveTerm) {
      powerTerm = x*x*powerTerm; 
      factorialTerm = i*(i-1)*factorialTerm;    
      result += (powerTerm / factorialTerm);
      positiveTerm = 0;
    } else {
      powerTerm = x*x*powerTerm; 
      factorialTerm = i*(i-1)*factorialTerm;   
      result -= (powerTerm / factorialTerm);
      positiveTerm = 1;
    }

  }

  return result;
}

float normalizeRad(float x) {
  /* 
    Taking input a float and normalizes it to range [-pi, pi]
  */

  if (x > pi) {
    while (x > pi) {
      x -= 2 * pi;   // If we have positive float greater than pi, keep subtracting until it yes than pi
    }
  } else if (x < -pi) {
    while (x < -pi) {
      x += 2 * pi;
    }
  }

  return x;
}