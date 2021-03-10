#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>

#include "timer.h"

namespace NM {

auto f = [](double x) { return std::pow(x, 4) - 2 * x + 1; };

class Integral {
 public:
  Integral(std::function<double(double)> f) : m_f(f) {}
  ~Integral() = default;

  virtual double integrate(double a, double b, size_t nSteps) const = 0;

 protected:
  std::function<double(double)> m_f;
};

class Trapezium : public Integral {
 public:
  Trapezium(std::function<double(double)> f) : Integral(f) {}
  ~Trapezium() = default;

  double integrate(double a, double b, size_t nSteps) const {
    const double h = (b - a) / (double)(nSteps);
    double value = 0.5 * m_f(a) + 0.5 * m_f(b);
    for (size_t k = 1; k < nSteps; ++k) {
      value += m_f(a + k * h);
    }
    return h * value;
  }
};

}  // namespace NM

int main() {
  NM::Timer timer;
  const double accuracyRequired = 1e-30;
  double currentAccuracy = 1e10;
  double currentStepOrder = 5;
  NM::Trapezium trapezium(NM::f);
  while (currentAccuracy > accuracyRequired && currentStepOrder < 30) {
    const double nSteps = std::pow(2, currentStepOrder);
    const double I = trapezium.integrate(0, 2, nSteps);
    currentAccuracy = std::abs(I - 4.4);
    std::cout << currentStepOrder << "\t" << std::setprecision(5) << I << "\t" << currentAccuracy << "\n";
    currentStepOrder++;
  }
  return 0;
}