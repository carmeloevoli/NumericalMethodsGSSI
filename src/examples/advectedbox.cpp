#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include "gsl.h"
#include "timer.h"
#include "utils.h"

namespace NM {

constexpr double u = 0.1;
constexpr double H = 1.;

std::vector<double> Advection(std::vector<double> z, std::vector<double> f, double dt) {
  auto N = f.size();
  std::vector<double> fNew(N);
  std::fill(fNew.begin(), fNew.end(), 0.);
  const auto dz = z[1] - z[0];

  std::vector<double> centralDiagonal(N - 2);
  std::vector<double> rhs(N - 2);
  std::vector<double> solution(N - 2);
  std::vector<double> lowerDiagonal(N - 3);
  std::vector<double> upperDiagonal(N - 3);

  const auto gamma = u * dt / dz;

  for (size_t i = 1; i < N - 1; ++i) {
    centralDiagonal.at(i - 1) = 1. + 0.5 * gamma;
    if (i != 1) lowerDiagonal.at(i - 2) = -0.5 * gamma;
    if (i != N - 2) upperDiagonal.at(i - 1) = 0.;
  }

  for (size_t i = 1; i < N - 1; ++i) {
    const auto Q = 0.;
    rhs.at(i - 1) = dt * Q;
    rhs.at(i - 1) += 0.5 * gamma * f.at(i - 1) + (1. - 0.5 * gamma) * f.at(i);
  }

  GSL::gsl_linalg_solve_tridiag(centralDiagonal, upperDiagonal, lowerDiagonal, rhs, solution);

  for (size_t i = 1; i < N - 1; ++i) fNew[i] = solution[i - 1];
  return fNew;
}

std::vector<double> AdvectionLax(std::vector<double> z, std::vector<double> f, double dt) {
  auto N = f.size();
  std::vector<double> fNew(N);
  std::fill(fNew.begin(), fNew.end(), 0.);
  const auto dz = z[1] - z[0];

  std::vector<double> centralDiagonal(N - 2);
  std::vector<double> rhs(N - 2);
  std::vector<double> solution(N - 2);
  std::vector<double> lowerDiagonal(N - 3);
  std::vector<double> upperDiagonal(N - 3);

  const auto gamma = u * dt / 4. / dz;

  for (size_t i = 1; i < N - 1; ++i) {
    centralDiagonal.at(i - 1) = 1.;
    if (i != 1) lowerDiagonal.at(i - 2) = -gamma;
    if (i != N - 2) upperDiagonal.at(i - 1) = gamma;
  }

  for (size_t i = 1; i < N - 1; ++i) {
    const auto Q = 0.;
    rhs.at(i - 1) = dt * Q;
    rhs.at(i - 1) += gamma * f.at(i - 1) + f.at(i) - gamma * f.at(i + 1);
  }

  GSL::gsl_linalg_solve_tridiag(centralDiagonal, upperDiagonal, lowerDiagonal, rhs, solution);

  for (size_t i = 1; i < N - 1; ++i) fNew[i] = solution[i - 1];
  return fNew;
}

std::vector<double> initialCondition(std::vector<double> z) {
  std::vector<double> f;
  for (auto i = 0; i < z.size(); ++i) {
    if (z[i] >= 0.1 && z[i] <= 0.2)
      f.push_back(1.);
    else
      f.push_back(0.);
  }
  return f;
}

}  // namespace NM

double runSim(int zSteps, double dt, std::string outputFilename) {
  NM::Timer timer(outputFilename);
  const auto z = NM::linspace(0., NM::H, zSteps);
  const auto dz = z[1] - z[0];

  auto f = NM::initialCondition(z);

  std::cout << "Courant : " << NM::u * dt / dz << "\n";
  for (size_t i = 0; i < 10; ++i) {
    std::cout << i << "\n";
    NM::saveSolution(z, f, i, outputFilename);
    for (size_t j = 0; j < (size_t)(1. / dt); ++j) {
      f = NM::AdvectionLax(z, f, dt);
    }
  }

  return 0;
}

int main() {
  std::cout << "t_adv : " << NM::H / NM::u << "\n";
  runSim(1001, 0.01, "box_3");
  runSim(10001, 0.001, "box_4");
  runSim(100001, 0.0001, "box_5");
  return 0;
}