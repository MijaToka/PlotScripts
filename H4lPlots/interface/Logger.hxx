#ifndef LOGGER_HXX
#define LOGGER_HXX

#include <iostream>

namespace logging {
inline bool verbose = false;
}

#define LOG(msg)                                                               \
  do {                                                                         \
    if (logging::verbose)                                                      \
      std::cout << "[INFO]\t" << msg << std::endl;                             \
  } while (0)

#define LOGERR(msg)                                                            \
  do {                                                                         \
    if (logging::verbose)                                                      \
      std::cerr << "[ERROR]\t" << msg << std::endl;                            \
  } while (0)

#endif
