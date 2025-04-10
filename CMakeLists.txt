##
#  CMake script for the step-1 tutorial program:
##

# Set the name of the project and target:
set(TARGET "main")

FILE(GLOB_RECURSE TARGET_SRC "source/*.cpp")
INCLUDE_DIRECTORIES(include)

# set(TARGET_SRC
#   ${TARGET}.cpp
#   )

# Usually, you will not need to modify anything beyond this point...

cmake_minimum_required(VERSION 3.13.4)

find_package(deal.II 9.6.0
  HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
  )
if(NOT ${deal.II_FOUND})
  message(FATAL_ERROR "\n"
    "*** Could not locate a (sufficiently recent) version of deal.II. ***\n\n"
    "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake\n"
    "or set an environment variable \"DEAL_II_DIR\" that contains this path."
    )
endif()

deal_ii_initialize_cached_variables()
project(${TARGET})
deal_ii_invoke_autopilot()

include(CMakePrintHelpers)
find_package(ideal.II 0.9.0 
  HINTS ${ideal.II_DIR} ${IDEAL_II_DIR} ../ ../../ $ENV{IDEAL_II_DIR}
)
 
if(NOT ${ideal.II_FOUND})
  message(FATAL_ERROR "\n"
    "*** Could not locate a (sufficiently recent) version of ideal.II. ***\n\n"
    "You may want to either pass a flag -DIDEAL_II_DIR=/path/to/ideal.II to cmake\n"
    "or set an environment variable \"IDEAL_II_DIR\" that contains this path."
  )
endif()

target_link_libraries(${TARGET} ideal.II::ideal.II)