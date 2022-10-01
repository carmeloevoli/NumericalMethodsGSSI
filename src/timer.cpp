#include "timer.h"

namespace NM {

Timer::Timer(std::string message) : m_message(message), m_start(AwesomeClock::now()) {}

Timer::~Timer() {
  m_end = std::chrono::high_resolution_clock::now();
  m_duration = m_end - m_start;
  std::cerr << m_message << " " << m_duration.count() << " s.\n";
}

}  // namespace NM
