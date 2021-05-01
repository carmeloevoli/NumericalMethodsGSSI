#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include "utils.h"

namespace NM {

double wave(double x) { return std::sin(8.0 * x * M_PI); }

std::vector<double> computeInitialCondition(std::vector<double> x) {
  const auto xSize = x.size();
  std::vector<double> ic;
  ic.reserve(xSize);
  for (auto xi : x) {
    const auto value = wave(xi);
    ic.emplace_back(value);
  }
  return ic;
}

std::vector<double> LaxFriedrichs(std::vector<double> f, double alpha) {
  auto N = f.size();
  std::vector<double> fNew(N);
  double fNext = 0;
  for (size_t i = 0; i < N; ++i) {
    if (i == 0) {
      fNext = 0.5 * (f[1] + f[N - 1]) - 0.5 * alpha * (f[1] - f[N - 1]);
    } else if (i == N - 1) {
      fNext = 0.5 * (f[0] + f[N - 2]) - 0.5 * alpha * (f[0] - f[N - 2]);
    } else {
      fNext = 0.5 * (f[i + 1] + f[i - 1]) - 0.5 * alpha * (f[i + 1] - f[i - 1]);
    }
    fNew.at(i) = fNext;
  }
  return fNew;
}

std::vector<double> Upwind(std::vector<double> f, double alpha) {
  auto N = f.size();
  std::vector<double> fNew(N);
  double fNext = 0;
  for (size_t i = 0; i < N; ++i) {
    if (i == 0) {
      fNext = f[0] - alpha * (f[0] - f[N - 1]);
    } else {
      fNext = f[i] - alpha * (f[i] - f[i - 1]);
    }
    fNew.at(i) = fNext;
  }
  return fNew;
}

std::vector<double> LaxWendroff(std::vector<double> f, double alpha) {
  auto N = f.size();
  std::vector<double> fNew(N);
  double fNext = 0;
  for (size_t i = 0; i < N; ++i) {
    if (i == 0) {
      const double fHalfDown = 0.5 * (f[0] + f[N - 1]) - 0.5 * alpha * (f[0] - f[N - 1]);
      const double fHalfUp = 0.5 * (f[1] + f[0]) - 0.5 * alpha * (f[1] - f[0]);
      fNext = f[0] - alpha * (fHalfUp - fHalfDown);
    } else if (i == N - 1) {
      const double fHalfDown = 0.5 * (f[N - 1] + f[N - 2]) - 0.5 * alpha * (f[N - 1] - f[N - 2]);
      const double fHalfUp = 0.5 * (f[0] + f[N - 1]) - 0.5 * alpha * (f[0] - f[N - 1]);
      fNext = f[N - 1] - alpha * (fHalfUp - fHalfDown);
    } else {
      const double fHalfDown = 0.5 * (f[i] + f[i - 1]) - 0.5 * alpha * (f[i] - f[i - 1]);
      const double fHalfUp = 0.5 * (f[i + 1] + f[i]) - 0.5 * alpha * (f[i + 1] - f[i]);
      fNext = f[i] - alpha * (fHalfUp - fHalfDown);
    }
    fNew.at(i) = fNext;
  }
  return fNew;
}

double computeRMSError(std::vector<double> x, std::vector<double> f, double v) {
  auto N = x.size();
  double error = 0;
  for (size_t i = 0; i < N; ++i) error += std::pow(f[i] - wave(x[i]), 2.0);
  auto rmsError = std::sqrt(error / (double)N);
  return rmsError;
}

double computeAverageError(std::vector<double> x, std::vector<double> f, double v) {
  auto N = x.size();
  double error = 0;
  for (size_t i = 1; i < N - 1; ++i) error += std::abs(f[i] - wave(x[i]));
  auto meanError = error / (double)N;
  return meanError;
}

}  // namespace NM

enum Algorithms { UPWIND, LAXF, LAXW };

double runSim(Algorithms algorithm, double dt, double tMax, int xOrder, std::string outputFilename) {
  const size_t xSteps = std::pow(2, xOrder);
  const double xMin = 0;
  const double xMax = 1;
  const auto x = NM::linspace(xMin, xMax, xSteps);
  const auto dx = x[1] - x[0];

  const double v = 1.;
  const double alpha = dt * v / dx;

  if (alpha > 0.99) std::cout << "!CFL violated with alpha = " << alpha << "\n";

  const auto initialCondition = NM::computeInitialCondition(x);

  std::vector<double> f;
  std::copy(initialCondition.begin(), initialCondition.end(), std::back_inserter(f));

  NM::saveSolution(x, f, 0, outputFilename);

  double t = tMax;

  if (algorithm == UPWIND) {
    while (t > dt) {
      f = NM::Upwind(f, alpha);
      t -= dt;
    }
    f = NM::Upwind(f, v * t / dx);
  } else if (algorithm == LAXF) {
    while (t > dt) {
      f = NM::LaxFriedrichs(f, alpha);
      t -= dt;
    }
    f = NM::LaxFriedrichs(f, v * t / dx);
  } else if (algorithm == LAXW) {
    while (t > dt) {
      f = NM::LaxWendroff(f, alpha);
      t -= dt;
    }
    f = NM::LaxWendroff(f, v * t / dx);
  } else {
    std::cout << "No algorithm found!\n";
  }

  NM::saveSolution(x, f, 1, outputFilename);

  return NM::computeAverageError(x, f, v);
}

int main() {
  {
    runSim(UPWIND, 0.001, 1, 7, "upwind_1_7");
    runSim(UPWIND, 0.001, 1, 8, "upwind_1_8");
    runSim(UPWIND, 0.001, 1, 9, "upwind_1_9");
    runSim(LAXF, 0.001, 1, 7, "laxf_1_7");
    runSim(LAXF, 0.001, 1, 8, "laxf_1_8");
    runSim(LAXF, 0.001, 1, 9, "laxf_1_9");
    runSim(LAXW, 0.001, 1, 7, "laxw_1_7");
    runSim(LAXW, 0.001, 1, 8, "laxw_1_8");
    runSim(LAXW, 0.001, 1, 9, "laxw_1_9");
  }

  {
    for (size_t i = 20; i > 4; --i) {
      std::cout << std::fixed << i << " ";
      std::cout << std::scientific;
      std::cout << runSim(UPWIND, .5e-6, 1, i, "test") << " ";
      std::cout << runSim(LAXF, .5e-6, 1, i, "test") << " ";
      std::cout << runSim(LAXW, .5e-6, 1, i, "test") << " ";
      std::cout << "\n";
    }
  }

  return 0;
}