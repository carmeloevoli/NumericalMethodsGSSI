#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <memory>

#include "timer.h"

namespace NM {

struct f {
  double operator()(double x) { return std::pow(x, 4) - 2 * x + 1; }
};

template <class Function>
double integrate(Function f, double a, double b, size_t nSteps) {
  const double h = (b - a) / (double)(nSteps);
  double value = 0.5 * f(a) + 0.5 * f(b);
  for (size_t k = 1; k < nSteps; ++k) {
    value += f(a + k * h);
  }
  return h * value;
}

}  // namespace NM

int main() {
  NM::Timer timer;
  const double accuracyRequired = 1e-30;
  double currentAccuracy = 1e10;
  double currentStepOrder = 5;
  while (currentAccuracy > accuracyRequired && currentStepOrder < 30) {
    const double nSteps = std::pow(2, currentStepOrder);
    const double I = integrate(NM::f(), 0.0, 2.0, nSteps);
    currentAccuracy = std::abs(I - 4.4);
    std::cout << currentStepOrder << "\t" << std::setprecision(5) << I << "\t" << currentAccuracy << "\n";
    currentStepOrder++;
  }
  return 0;
}