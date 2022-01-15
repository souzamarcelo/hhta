#pragma once

#include <random>
#include <chrono>
#include <iostream>
#include <iomanip>

extern bool verbose;
extern bool record;
extern unsigned timelimit;

extern std::discrete_distribution<> rank_proportional; // rank-proportional distribution over variables

inline bool termination(int bestValue, int worstValue, unsigned steps, int target, std::chrono::system_clock::time_point start, std::chrono::system_clock::time_point& last_report) {
  if (std::chrono::system_clock::now()-last_report>std::chrono::duration<int>(3)) {
    last_report = std::chrono::system_clock::now();
    if (verbose)
      std::cerr << "RC" << std::setw(7) << steps << " " << std::setw(3) << std::chrono::duration_cast<std::chrono::seconds>(last_report-start).count() << " " << std::setw(7) << bestValue << " " << std::setw(7) << worstValue << std::endl;
  }
  
  if (bestValue <= target)
    return true;
  
  if (std::chrono::system_clock::now()-start>std::chrono::duration<int>(timelimit))
    return true;
  return false;
}

