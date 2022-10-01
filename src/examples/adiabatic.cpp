#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include "gsl.h"
#include "timer.h"
#include "utils.h"

namespace NM {

constexpr double Q_0 = 1.;
constexpr double alpha = 4.;
constexpr double beta = 10.;  // s^-1
constexpr double p_0 = 1.;    // GV

double f(double p) {
  constexpr double C = 3. * Q_0 / alpha / beta;
  return C * (1. - std::pow(p / p_0, -alpha));
};

std::vector<double> Upwind(std::vector<double> p, std::vector<double> f, double dt) {
  auto N = f.size();
  std::vector<double> fNew(N);
  double fNext = 0;
  const auto dlgp = p[1] / p[0];
  for (size_t i = 1; i < N; ++i) {
    const auto Q = Q_0 * 0.5 * (std::pow(p.at(i) / p_0, -alpha) + std::pow(p.at(i - 1) / p_0, -alpha));
    const auto C = beta / 3. * (dlgp / (dlgp - 1.));
    fNext = f.at(i) + dt * (Q - C * (f.at(i) - f.at(i - 1)));
    fNew.at(i) = fNext;
  }
  return fNew;
}

std::vector<double> CrankNicolson(std::vector<double> p, std::vector<double> f, double dt) {
  auto N = f.size();
  std::vector<double> fNew(N);
  const auto dlgp = p[1] / p[0];

  std::vector<double> centralDiagonal(N - 1);
  std::vector<double> rhs(N - 1);
  std::vector<double> lowerDiagonal(N - 2);
  std::vector<double> upperDiagonal(N - 2);

  const auto C = beta / 3. * (dlgp / (dlgp - 1.));
  std::fill(lowerDiagonal.begin(), lowerDiagonal.end(), -dt / 2. * C);
  std::fill(centralDiagonal.begin(), centralDiagonal.end(), 1. + dt / 2. * C);
  std::fill(upperDiagonal.begin(), upperDiagonal.end(), 0.);

  for (size_t i = 1; i < N; ++i) {
    const auto Q = Q_0 * 0.5 * (std::pow(p.at(i) / p_0, -alpha) + std::pow(p.at(i - 1) / p_0, -alpha));
    rhs[i - 1] = dt / 2. * C * f[i - 1] + (1. - dt / 2. * C) * f[i] + dt * Q;
  }

  std::vector<double> solution(N - 1);
  GSL::gsl_linalg_solve_tridiag(centralDiagonal, upperDiagonal, lowerDiagonal, rhs, solution);

  for (size_t i = 1; i < N; ++i) fNew[i] = solution[i - 1];
  fNew[0] = 0.;
  return fNew;
}

std::vector<double> CrankNicolsonLog(std::vector<double> y, std::vector<double> f, double dt) {
  auto N = f.size();
  std::vector<double> fNew(N);
  const auto dy = y[1] - y[0];

  std::vector<double> centralDiagonal(N - 1);
  std::vector<double> rhs(N - 1);
  std::vector<double> lowerDiagonal(N - 2);
  std::vector<double> upperDiagonal(N - 2);

  const auto C = beta / 3. / dy;
  std::fill(lowerDiagonal.begin(), lowerDiagonal.end(), -dt / 2. * C);
  std::fill(centralDiagonal.begin(), centralDiagonal.end(), 1. + dt / 2. * C);
  std::fill(upperDiagonal.begin(), upperDiagonal.end(), 0.);

  for (size_t i = 1; i < N; ++i) {
    const auto p = std::exp(y.at(i));
    const auto pLow = std::exp(y.at(i - 1));
    const auto Q = Q_0 * 0.5 * (std::pow(p / p_0, -alpha) + std::pow(pLow / p_0, -alpha));
    rhs[i - 1] = dt / 2. * C * f[i - 1] + (1. - dt / 2. * C) * f[i] + dt * Q;
  }

  std::vector<double> solution(N - 1);
  GSL::gsl_linalg_solve_tridiag(centralDiagonal, upperDiagonal, lowerDiagonal, rhs, solution);

  for (size_t i = 1; i < N; ++i) fNew[i] = solution[i - 1];
  fNew[0] = 0.;
  return fNew;
}

}  // namespace NM

void printSolution(double pMin, double pMax, std::string outputFilename) {
  const auto p = NM::logspace(pMin, pMax, 10000);
  std::vector<double> f_p;
  for (auto& pi : p) f_p.push_back(NM::f(pi));
  NM::saveSolution(p, f_p, 0, outputFilename);
}

double runSim(double pMin, double pMax, int xOrder, std::string outputFilename, bool doCN = true) {
  printSolution(pMin, pMax, "solution");
  NM::Timer timer(outputFilename);
  const size_t pSteps = std::pow(2, xOrder);
  const auto p = NM::logspace(pMin, pMax, pSteps);
  std::vector<double> f(pSteps);
  std::fill(f.begin(), f.end(), 0.);

  const double dt = 0.001;
  for (size_t i = 0; i < 10; ++i) {
    if (doCN) {
      for (size_t j = 0; j < 500; ++j) f = NM::CrankNicolson(p, f, dt);
    } else {
      for (size_t j = 0; j < 500; ++j) f = NM::Upwind(p, f, dt);
    }
    NM::saveSolution(p, f, i, outputFilename);
  }

  return 0;
}

double runSimLog(double pMin, double pMax, int xOrder, std::string outputFilename) {
  printSolution(pMin, pMax, "solution");
  NM::Timer timer(outputFilename);
  const size_t pSteps = std::pow(2, xOrder);
  const auto y = NM::linspace(std::log(pMin), std::log(pMax), pSteps);
  std::vector<double> f(pSteps);
  std::fill(f.begin(), f.end(), 0.);

  const double dt = 0.001;
  for (size_t i = 0; i < 10; ++i) {
    for (size_t j = 0; j < 500; ++j) f = NM::CrankNicolsonLog(y, f, dt);
    NM::saveSolution(y, f, i, outputFilename);
  }

  return 0;
}

int main() {
  runSim(NM::p_0, 1e4 * NM::p_0, 7, "CN", true);
  runSim(NM::p_0, 1e4 * NM::p_0, 7, "upwind", false);
  runSimLog(NM::p_0, 1e4 * NM::p_0, 6, "CNLog");
  return 0;
}
