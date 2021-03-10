#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

namespace NM {

auto func = [](double x) { return std::pow(x, 4) - 2 * x + 1; };

struct c_unique {
  int current;
  c_unique() { current = -1; }
  int operator()() { return ++current; }
} UniqueNumber;

double integrate(double a, double b, size_t nSteps) {
  std::vector<double> fArray(nSteps);
  std::generate(fArray.begin(), fArray.end(), UniqueNumber);

  //  for (auto x_i : x) std::cout << x_i << "\n";
  return 0;
}

}  // namespace NM

int main() { NM::integrate(0, 2, 32); }