#ifndef INCLUDE_TIMER_H_
#define INCLUDE_TIMER_H_

#include <chrono>
#include <iostream>

namespace NM {

using AwesomeClock = std::chrono::high_resolution_clock;

class Timer {
 public:
  Timer();
  ~Timer();

 protected:
  std::chrono::time_point<AwesomeClock> m_start;
  std::chrono::time_point<AwesomeClock> m_end;
  std::chrono::duration<double> m_duration;
};

}  // namespace NM

#endif