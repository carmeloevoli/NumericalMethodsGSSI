#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include "gsl.h"
#include "timer.h"
#include "utils.h"

namespace NM {

constexpr double kpc = 1;
constexpr double Myr = 1;

constexpr double s = Myr / 3.154e13;
constexpr double cm = kpc / 3.086e21;

constexpr double Q_0 = 1.;
constexpr double diskSize = 0.1 * kpc;
constexpr double HaloSize = 10. * kpc;
constexpr double u = 4. * 1e6 * cm / s;
constexpr double D = 3e28 * cm * cm / s;
constexpr double beta = u * HaloSize / D;

// double f(double z) {
//   constexpr double C = diskSize * Q_0 / u;
//   return C * std::exp(beta) * (1. - std::exp(-beta * (1. - std::fabs(z) / HaloSize)));
// }

double f(double z) { return diskSize * HaloSize * Q_0 / D * (1. - std::fabs(z) / HaloSize); }

std::vector<double> CrankNicolson(std::vector<double> z, std::vector<double> f, double dt) {
  auto N = f.size();
  std::vector<double> fNew(N);
  std::fill(fNew.begin(), fNew.end(), 0.);
  const auto dz = z[1] - z[0];

  std::vector<double> centralDiagonal(N - 2);
  std::vector<double> rhs(N - 2);
  std::vector<double> lowerDiagonal(N - 3);
  std::vector<double> upperDiagonal(N - 3);

  const auto alpha = D * dt / std::pow(dz, 2.);

  for (size_t i = 1; i < N - 1; ++i) {
    centralDiagonal.at(i - 1) = 1. + alpha;
    if (i != 1) lowerDiagonal.at(i - 2) = -0.5 * alpha;
    if (i != N - 2) upperDiagonal.at(i - 1) = -0.5 * alpha;
  }

  for (size_t i = 1; i < N - 1; ++i) {
    const auto Q = 2. * diskSize * Q_0 * NM::Gaussian1D(std::fabs(z.at(i)), diskSize);
    rhs.at(i - 1) = dt * Q;
    rhs.at(i - 1) += 0.5 * alpha * f.at(i + 1) + (1. - alpha) * f.at(i) + 0.5 * alpha * f.at(i - 1);
  }

  std::vector<double> solution(N - 2);
  GSL::gsl_linalg_solve_tridiag(centralDiagonal, upperDiagonal, lowerDiagonal, rhs, solution);

  for (size_t i = 1; i < N - 1; ++i) fNew[i] = solution[i - 1];
  return fNew;
}

}  // namespace NM

void printSolution(std::string outputFilename) {
  const auto z = NM::linspace(-NM::HaloSize, NM::HaloSize, 10000);
  std::vector<double> f_z;
  for (auto& z_i : z) f_z.push_back(NM::f(z_i));
  NM::saveSolution(z, f_z, 0, outputFilename);
}

double runSim(int zOrder, std::string outputFilename) {
  printSolution("solution");
  NM::Timer timer(outputFilename);
  const size_t zSteps = std::pow(2, zOrder) + 1;
  const auto z = NM::linspace(-NM::HaloSize, NM::HaloSize, zSteps);

  std::vector<double> f(zSteps);
  std::fill(f.begin(), f.end(), 0.);

  const double dt = 0.01;
  for (size_t i = 0; i < 10; ++i) {
    for (size_t j = 0; j < 20000; ++j) f = NM::CrankNicolson(z, f, dt);
    NM::saveSolution(z, f, i, outputFilename);
  }

  return 0;
}

int main() {
  std::cout << "t_adv : " << NM::HaloSize / NM::u << "\n";
  std::cout << "t_dif : " << std::pow(NM::HaloSize, 2) / NM::D << "\n";
  std::cout << "beta : " << NM::beta << "\n";

  runSim(9, "galaxy");
  return 0;
}
