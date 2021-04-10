#include <algorithm>
#include <cmath>
#include <iomanip>
#include <random>
#include <vector>

#include "timer.h"

using point = std::pair<double, double>;

std::pair<double, double> randPoint() {
  static std::random_device rd;
  static std::mt19937 mt_eng(rd());  // mersenne-twister engine initialised with seed
  std::uniform_real_distribution<> dist(0., 1.);
  return {dist(mt_eng), dist(mt_eng)};
}

double f(double x) { return 0.5 * std::sin(x * 8.0 * M_PI) + 0.5; }
bool isBelowCurve(point p) { return f(p.first) < p.second; }

int main() {
  NM::Timer timer;
  int N = 10000;
  auto points = std::vector<point>();
  points.resize(N);
  std::generate(points.begin(), points.end(), randPoint);
  //  for (auto point : points) std::cout << point.first << " " << point.second << "\n";
  int goodPoints = std::count_if(points.begin(), points.end(), isBelowCurve);
  std::cout << goodPoints << " " << std::setprecision(3) << goodPoints / (double)N << "\n";
}