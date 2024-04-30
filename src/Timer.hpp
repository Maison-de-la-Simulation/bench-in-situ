#ifndef EULER_SRC_TIMER_HPP_
#define EULER_SRC_TIMER_HPP_

#include <chrono>

struct PerformanceTimer {
  PerformanceTimer() : time_spent_in_io(0), time_spent_in_compute(0) {
  }

  friend std::ostream& operator<<(std::ostream& os, const PerformanceTimer& obj) {
    os << "[RESULT] I/O: " << obj.time_spent_in_compute.count() << " s" << std::endl;
    os << "[RESULT] Compute: " << obj.time_spent_in_io.count() << " s" << std::endl;
    return os;
  }

  std::chrono::duration<double> time_spent_in_io;
  std::chrono::duration<double> time_spent_in_compute;
};



#endif //EULER_SRC_TIMER_HPP_
