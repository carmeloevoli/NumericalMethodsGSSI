#include "utils.h"

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