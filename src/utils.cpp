#include "utils.h"

#include <cmath>
#include <fstream>
#include <iostream>

namespace NM {

std::vector<double> linspace(double xMin, double xMax, size_t xSize) {
  std::vector<double> x;
  x.reserve(xSize);
  const auto h = (xMax - xMin) / (double)(xSize - 1);
  for (size_t i = 0; i < xSize; ++i) {
    const auto xi = xMin + (double)i * h;
    x.emplace_back(xi);
  }
  return x;
}

std::vector<double> logspace(double xMin, double xMax, size_t xSize) {
  std::vector<double> x;
  x.reserve(xSize);
  const auto delta_log = std::exp(std::log(xMax / xMin) / (xSize - 1));
  for (size_t i = 0; i < xSize; ++i) {
    const auto xi = std::exp(std::log(xMin) + (double)i * std::log(delta_log));
    x.emplace_back(xi);
  }
  return x;
}

double Gaussian1D(double z, double sigma) {
  return 1. / sigma / std::sqrt(2. * M_PI) * std::exp(-0.5 * std::pow(z / sigma, 2.));
}

void saveSolution(std::vector<double> x, std::vector<double> f, int counter, std::string filename) {
  std::ofstream txtfile;
  txtfile.open("output/" + filename + "_" + std::to_string(counter) + ".txt");
  if (txtfile.is_open()) {
    for (size_t i = 0; i < x.size(); ++i) txtfile << x.at(i) << " " << f.at(i) << "\n";
    txtfile.close();
  } else
    std::cout << "!Unable to open file";
}

}  // namespace NM