#include "timer.h"
#include "vector3.h"

namespace NM {

struct LorenzParams {
  double sigma;
  double r;
  double b;
};

class Func {
 public:
  Func(double sigma, double r, double b) : m_params({sigma, r, b}) {}

  Vector3d eval(Vector3d v) {
    const auto dxdt = m_params.sigma * (v.y - v.x);
    const auto dydt = m_params.r * v.x - v.y - v.x * v.z;
    const auto dzdt = v.x * v.y - m_params.b * v.z;
    return Vector3d(dxdt, dydt, dzdt);
  }

 protected:
  LorenzParams m_params;
};

class Lorenz {
 public:
  Lorenz(Func f) : m_f(f) {}
  void setInitialCondition(Vector3d initialCondition) { m_v = initialCondition; };
  void evolveInTime(double beginTime, double endTime, size_t N) {
    const auto h = (endTime - beginTime) / (double)N;
    for (size_t i = 0; i < N; ++i) {
      const auto v = m_v;
      const auto k_1 = h * m_f.eval(v);
      const auto k_2 = h * m_f.eval(v + 0.5 * k_1);
      const auto k_3 = h * m_f.eval(v + 0.5 * k_2);
      const auto k_4 = h * m_f.eval(v + 0.5 * k_3);
      m_v = v + (1. / 6.) * (k_1 + 2. * k_2 + 2. * k_3 + k_4);
      std::cout << beginTime + (double)i * h << " " << m_v << "\n";
    }
  }

 protected:
  Vector3d m_v;
  Func m_f;
};

}  // namespace NM

int main() {
  NM::Timer timer;
  NM::Func func(10., 28., 8. / 3.);
  NM::Lorenz lorenz(func);
  lorenz.setInitialCondition(NM::Vector3d(0., 1., 0.));
  lorenz.evolveInTime(0., 50., 100000);
  return 0;
}