#ifndef INCLUDE_TIMER_H_
#define INCLUDE_TIMER_H_

#include <chrono>
#include <iostream>

namespace NM {

using AwesomeClock = std::chrono::high_resolution_clock;

class Timer {
 protected:
  std::string m_message;
  std::chrono::time_point<AwesomeClock> m_start;
  std::chrono::time_point<AwesomeClock> m_end;
  std::chrono::duration<double> m_duration;

 public:
  explicit Timer(std::string message);
  ~Timer();
};

}  // namespace NM

#endif