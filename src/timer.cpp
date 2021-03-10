#include "timer.h"

namespace NM {
Timer::Timer() { m_start = AwesomeClock::now(); }

Timer::~Timer() {
  m_end = std::chrono::high_resolution_clock::now();
  m_duration = m_end - m_start;
  std::cerr << "Timer took : " << m_duration.count() << " s." << std::endl;
}

}  // namespace NM
