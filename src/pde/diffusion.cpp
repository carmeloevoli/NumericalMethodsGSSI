#include <cmath>
#include <vector>

#include "gsl.h"
#include "utils.h"

#define pow2(A) ((A) * (A))

namespace NM {

auto fa = [](double x) { return std::cos(0.5 * M_PI * x); };

std::vector<double> sourceTerm(std::vector<double> x, double D) {
  std::vector<double> Q;
  Q.reserve(x.size());
  for (auto xi : x) Q.emplace_back(0.25 * D * pow2(M_PI) * std::cos(0.5 * M_PI * xi));
  return Q;
}

std::vector<double> CrankNicolson(std::vector<double> f, std::vector<double> Q, double alpha, double dt) {
  auto N = f.size();
  std::vector<double> fNew(N);

  std::vector<double> centralDiagonal(N - 2);
  std::vector<double> rhs(N - 2);
  std::vector<double> lowerDiagonal(N - 3);
  std::vector<double> upperDiagonal(N - 3);

  std::fill(centralDiagonal.begin(), centralDiagonal.end(), 1. + alpha);
  std::fill(lowerDiagonal.begin(), lowerDiagonal.end(), -0.5 * alpha);
  std::fill(upperDiagonal.begin(), upperDiagonal.end(), -0.5 * alpha);

  for (size_t i = 1; i < N - 1; ++i) {
    rhs[i - 1] = dt * Q[i] + 0.5 * alpha * f[i + 1] + (1. - alpha) * f[i] + 0.5 * alpha * f[i - 1];
  }

  std::vector<double> solution(N - 2);
  GSL::gsl_linalg_solve_tridiag(centralDiagonal, upperDiagonal, lowerDiagonal, rhs, solution);

  for (size_t i = 1; i < N - 1; ++i) fNew[i] = solution[i - 1];
  return fNew;
}

double computeRMSError(std::vector<double> x, std::vector<double> f) {
  auto N = x.size();
  double error = 0;
  for (size_t i = 0; i < N; ++i) error += std::pow(f[i] - fa(x[i]), 2.0);
  auto rmsError = std::sqrt(error / (double)N);
  return rmsError;
}

}  // namespace NM

double runSimCN(double dt, double tMax, int xOrder, std::string outputFilename) {
  const size_t xSteps = std::pow(2, xOrder);
  const double xMin = -1;
  const double xMax = 1;
  const auto x = NM::linspace(xMin, xMax, xSteps);
  const auto dx = x[1] - x[0];

  const double D = 0.1;
  const double alpha = dt * D / pow2(dx);
  // if (alpha > 0.99) std::cout << "!CFL violated with alpha = " << alpha << "\n";

  auto Q = NM::sourceTerm(x, D);

  std::vector<double> f(xSteps);

  NM::saveSolution(x, f, 0, outputFilename);

  double t = tMax;

  while (t > dt) {
    f = NM::CrankNicolson(f, Q, alpha, dt);
    t -= dt;
  }
  // f = NM::CrankNicolson(f, v * t / dx);

  NM::saveSolution(x, f, 1, outputFilename);

  return NM::computeRMSError(x, f);
}

int main() {
  for (double t = 10; t < 100; t += 1) {
    std::cout << t << " ";
    std::cout << runSimCN(0.001, t, 8, "diffusion_test") << " ";
    std::cout << runSimCN(0.001, t, 9, "diffusion_test") << " ";
    std::cout << runSimCN(0.001, t, 10, "diffusion_test") << " ";
    std::cout << "\n";
  }

  //   std::cout <<
  //   std::cout << runSimCN(0.001, 20, 8, "diffusion_20_8") << "\n";
  //   std::cout << runSimCN(0.001, 50, 8, "diffusion_5_8") << "\n";
  //   std::cout << runSimCN(0.001, 100, 8, "diffusion_100_8") << "\n";
}