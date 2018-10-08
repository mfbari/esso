#ifndef _IZ_TIMER_HPP_
#define _IZ_TIMER_HPP_

#include <chrono>
#include <utility>
#include <iostream>

class iz_timer {
    std::chrono::high_resolution_clock::time_point last_time_point;
    std::chrono::duration<double> time_duration;
  public:
    iz_timer() : 
      last_time_point {std::chrono::high_resolution_clock::now()},
      time_duration {std::chrono::duration<double>::zero()}
      {}
    double time() {
      auto n = std::chrono::high_resolution_clock::now();
      time_duration = n - last_time_point;
      return std::chrono::duration_cast<std::chrono::milliseconds>(
          time_duration).count();
    }
    void reset() {
      last_time_point = std::chrono::high_resolution_clock::now();
      time_duration = std::chrono::duration<double>::zero();
    }
    friend std::ostream& operator<<(std::ostream& out, iz_timer t) {
      out << t.time();
      return out;
    }
};

#endif
