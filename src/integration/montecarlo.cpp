#include <algorithm>
#include <cmath>
#include <iomanip>
#include <random>
#include <vector>

#include "timer.h"

namespace NM {

using point = std::pair<double, double>;

std::pair<double, double> randPoint() {
  static std::random_device rd;
  static std::mt19937 mt_eng(rd());  // mersenne-twister engine initialised with seed
  std::uniform_real_distribution<> dist(0., 1.);
  return {dist(mt_eng), dist(mt_eng)};
}

double f(double x) { return 0.5 * std::sin(x * 8.0 * M_PI) + 0.5; }
bool isBelowCurve(point p) { return f(p.first) < p.second; }

double integrate(const int N) {
  auto points = std::vector<point>();
  points.resize(N);
  std::generate(points.begin(), points.end(), randPoint);
  int goodPoints = std::count_if(points.begin(), points.end(), isBelowCurve);
  return (double)goodPoints / (double)N;
}

}  // namespace NM

int main() {
  NM::Timer timer;
  int currentStepOrder = 5;
  while (currentStepOrder < 25) {
    const auto nSteps = std::pow(2, currentStepOrder);
    for (size_t i = 0; i < 50; ++i) {
      const double I = NM::integrate(nSteps);
      std::cout << currentStepOrder << "\t" << std::setprecision(10) << I << "\n";
    }
    currentStepOrder++;
  }
  return 0;
}