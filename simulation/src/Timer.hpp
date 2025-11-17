#ifndef EULER_SRC_TIMER_HPP_
#define EULER_SRC_TIMER_HPP_

#include <chrono>

struct PerformanceTimer {
  PerformanceTimer() : time_spent_in_io(0), time_spent_in_compute(0) {
  }

  friend std::ostream& operator<<(std::ostream& os, const PerformanceTimer& obj) {
    os << "[RESULT] I/O_time " << obj.time_spent_in_io.count() << " s" << std::endl;
    os << "[RESULT] Compute_time " << obj.time_spent_in_compute.count() << " s" << std::endl;
    os << "[RESULT] Compute_/_I/O_time " << (obj.time_spent_in_compute.count())/(obj.time_spent_in_io.count()) << " Ratio" << std::endl;
    return os;
  }

  std::chrono::duration<double> time_spent_in_io;
  std::chrono::duration<double> time_spent_in_compute;
};

struct DebugTimer {
  DebugTimer() : time_spent_in_write_before_checkpoint(0), time_spent_in_write_after_checkpoint(0), time_spent_in_write_after_xml(0) {
  }

  friend std::ostream& operator<<(std::ostream& os, const DebugTimer& obj) {
    os << "[RESULT] Write_before_checkpoint " << obj.time_spent_in_write_before_checkpoint.count() << " s" << std::endl;
    os << "[RESULT] Write_after_checkpoint " << obj.time_spent_in_write_after_checkpoint.count() << " s" << std::endl;
    os << "[RESULT] Write_after_xml " << obj.time_spent_in_write_after_xml.count() << " s" << std::endl;
    return os;
  }

  std::chrono::duration<double> time_spent_in_write_before_checkpoint;
  std::chrono::duration<double> time_spent_in_write_after_checkpoint;
  std::chrono::duration<double> time_spent_in_write_after_xml;
};



#endif //EULER_SRC_TIMER_HPP_
