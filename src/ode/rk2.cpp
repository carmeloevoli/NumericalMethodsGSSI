#include <cmath>
#include <iomanip>
#include <iostream>

#include "timer.h"

namespace NM {

double dydx(double x, double y) { return x / y; }
double asolution(double x) { return std::sqrt(x * x + 1); }

double rungeKutta2nd(double x0, double y0, double xMax, size_t N) {
  const auto h = (xMax - x0) / (double)N;
  auto x = x0;
  auto y = y0;
  double maxRelativeError = 0;
  for (size_t i = 0; i < N; ++i) {
    const auto k1 = h * dydx(x, y);
    const auto k2 = h * dydx(x + 0.5 * h, y + 0.5 * k1);
    y += k2;
    x += h;
    const auto relativeError = std::abs(y - asolution(x)) / asolution(x);
    maxRelativeError = std::max(maxRelativeError, relativeError);
  }
  return maxRelativeError;
}

}  // namespace NM

int main() {
  NM::Timer timer;
  const double x0 = 0, y0 = 1, xMax = 10;
  double currentStepOrder = 5;
  while (currentStepOrder < 30) {
    const double nSteps = std::pow(2, currentStepOrder);
    const double maxError = NM::rungeKutta2nd(x0, y0, xMax, nSteps);
    std::cout << currentStepOrder << "\t" << std::setprecision(5) << maxError << "\n";
    currentStepOrder++;
  }
  return 0;
}