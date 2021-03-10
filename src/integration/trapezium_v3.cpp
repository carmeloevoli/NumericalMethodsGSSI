#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <vector>

#include "timer.h"

namespace NM {

auto func = [](double x) { return std::pow(x, 4) - 2 * x + 1; };

double integrate(double a, double b, size_t nSteps) {
  std::vector<double> fArray(nSteps - 1);
  const auto h = (b - a) / (double)nSteps;
  std::size_t k(1);
  std::generate(fArray.begin(), fArray.end(), [&] { return func(a + h * (k++)); });
  return h * std::accumulate(fArray.begin(), fArray.end(), 0.5 * func(a) + 0.5 * func(b));
}

}  // namespace NM

int main() {
  NM::Timer timer;
  const double accuracyRequired = 1e-30;
  double currentAccuracy = 1e10;
  double currentStepOrder = 5;
  while (currentAccuracy > accuracyRequired && currentStepOrder < 30) {
    const double nSteps = std::pow(2, currentStepOrder);
    const double I = NM::integrate(0.0, 2.0, nSteps);
    currentAccuracy = std::abs(I - 4.4);
    std::cout << currentStepOrder << "\t" << std::setprecision(5) << I << "\t" << currentAccuracy << "\n";
    currentStepOrder++;
  }
  return 0;
}