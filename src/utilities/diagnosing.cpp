#include <cmath>
#include <iostream>
#include <limits>

template <typename T>
int machep(bool doPrint = false) {
  const T one = T(1);
  T epsilon = T(1);
  while (!(one == T(1.) + epsilon)) {
    if (doPrint) std::cout << std::log2(epsilon) << "\n";
    epsilon /= T(2.);
  }
  return (int)std::log2(epsilon);
}

template <typename T>
int maxexp(bool doPrint = false) {
  T value = T(2);
  int counter = 1;
  while (value != std::numeric_limits<T>::infinity()) {
    if (doPrint) std::cout << std::log2(value) << "\n";
    value *= T(2);
    counter++;
  }
  return counter - 1;
}

int main() {
  std::cout << machep<float>() << "\n";
  std::cout << machep<double>() << "\n";
  std::cout << machep<long double>() << "\n";
  std::cout << "\n";
  std::cout << maxexp<float>() << "\n";
  std::cout << maxexp<double>() << "\n";
  std::cout << maxexp<long double>() << "\n";
}
