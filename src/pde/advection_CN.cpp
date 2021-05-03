#include <cmath>
#include <vector>

#include "gsl.h"
#include "utils.h"

namespace NM {

double gaussian(double x, double mu, double sigma) { return std::exp(-0.5 * std::pow((x - mu) / sigma, 2.0)); }

std::vector<double> CrankNicolson(std::vector<double> f, double alpha) {
  auto N = f.size();
  std::vector<double> fNew(N);

  std::vector<double> centralDiagonal(N - 2);
  std::vector<double> rhs(N - 2);
  std::vector<double> lowerDiagonal(N - 3);
  std::vector<double> upperDiagonal(N - 3);

  std::fill(centralDiagonal.begin(), centralDiagonal.end(), 1.);
  std::fill(lowerDiagonal.begin(), lowerDiagonal.end(), -0.25 * alpha);
  std::fill(upperDiagonal.begin(), upperDiagonal.end(), +0.25 * alpha);

  for (size_t i = 1; i < N - 1; ++i) {
    rhs[i - 1] = f[i] - 0.25 * alpha * (f[i + 1] - f[i - 1]);
  }

  std::vector<double> solution(N - 2);
  GSL::gsl_linalg_solve_tridiag(centralDiagonal, upperDiagonal, lowerDiagonal, rhs, solution);

  for (size_t i = 1; i < N - 1; ++i) fNew[i] = solution[i - 1];
  fNew[0] = 0.;
  fNew[N - 1] = 0.;
  return fNew;
}

double computeRMSError(std::vector<double> x, std::vector<double> f, double v) {
  auto N = x.size();
  double error = 0;
  for (size_t i = 0; i < N; ++i) error += std::pow(f[i] - gaussian(x[i], 0.1, 0.01), 2.0);
  auto rmsError = std::sqrt(error / (double)N);
  return rmsError;
}

double computeAverageError(std::vector<double> x, std::vector<double> f, double v) {
  auto N = x.size();
  double error = 0;
  for (size_t i = 1; i < N - 1; ++i) error += std::abs(f[i] - gaussian(x[i], 0.1, 0.01));
  auto meanError = error / (double)N;
  return meanError;
}

}  // namespace NM

double runSimCN(double dt, double tMax, int xOrder, std::string outputFilename) {
  const size_t xSteps = std::pow(2, xOrder);
  const double xMin = 0;
  const double xMax = 1;
  const auto x = NM::linspace(xMin, xMax, xSteps);
  const auto dx = x[1] - x[0];

  const double v = 1.;
  const double alpha = dt * v / dx;

  if (alpha > 0.99) std::cout << "!CFL violated with alpha = " << alpha << "\n";

  std::vector<double> f;
  f.reserve(xSteps);
  for (auto xi : x) {
    const auto value = NM::gaussian(xi, 0.1, 0.01);
    f.emplace_back(value);
  }

  NM::saveSolution(x, f, 0, outputFilename);

  double t = tMax;

  while (t > dt) {
    f = NM::CrankNicolson(f, alpha);
    t -= dt;
  }
  f = NM::CrankNicolson(f, v * t / dx);

  NM::saveSolution(x, f, 1, outputFilename);

  return NM::computeAverageError(x, f, v);
}

int main() {
  runSimCN(0.001, 0.6, 8, "CN_0.6_8");
  runSimCN(0.001, 0.6, 9, "CN_0.6_9");
  runSimCN(0.001, 0.6, 10, "CN_0.6_10");
  runSimCN(0.001, 0.6, 11, "CN_0.6_11");
}