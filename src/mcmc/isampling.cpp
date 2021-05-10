#include <algorithm>
#include <cmath>
#include <iomanip>
#include <random>
#include <vector>

#include "timer.h"

#define pow2(A) ((A) * (A))

namespace NM {

double randFromDistribution() {
  static std::random_device rd;
  static std::mt19937 mt_eng(rd());  // mersenne-twister engine initialised with seed
  std::uniform_real_distribution<> dist(0., 1.);
  auto r = dist(mt_eng);
  return pow2(r);
}

template <class Function>
double meanValue(Function f, const double a, const double b, const int N) {
  auto x = std::vector<double>(N);
  std::generate(x.begin(), x.end(), randFromDistribution);
  std::vector<double> y;
  y.reserve(N);
  for (auto xi : x) y.emplace_back(f(xi));
  return (b - a) / (double)N * std::accumulate(y.begin(), y.end(), 0.);
}

}  // namespace NM

int main() {
  NM::Timer timer;
  const auto nDraws = 1000;
  // auto f = [](double x) { return x * x + 1; };
  auto f = [](double x) { return 1. / (std::exp(x) + 1.); };
  const double integrationNorm = 2;
  for (size_t i = 0; i < 100; ++i) {
    const double I = NM::meanValue(f, 0, 1, nDraws);
    std::cout << i << "\t" << std::setprecision(4) << integrationNorm * I << "\n";
  }
  return 0;
}