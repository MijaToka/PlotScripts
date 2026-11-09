#ifndef LOGGER_HXX
#define LOGGER_HXX

#include <iostream>

namespace logging {
inline char verbose = 2;
}

#define LOG(msg, level)                                                        \
  do {                                                                         \
    if (logging::verbose >= level)                                             \
      std::cout << "[INFO]" << msg << std::endl;                               \
  } while (0)

#define LOGERR(msg, level)                                                     \
  do {                                                                         \
    if (logging::verbose >= level)                                             \
      std::cerr << "[ERROR]" << msg << std::endl;                              \
  } while (0)

#endif
