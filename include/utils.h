#include <string>
#include <vector>

namespace NM {

std::vector<double> linspace(double xMin, double xMax, size_t xSize);
std::vector<double> logspace(double xMin, double xMax, size_t xSize);
void saveSolution(std::vector<double> x, std::vector<double> f, int counter, std::string filename);

}  // namespace NM