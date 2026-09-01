cmake_minimum_required(VERSION 3.13)
project(PROAID VERSION 0.0.0.1 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_FLAGS "-g -O3")

include(GNUInstallDirs)
message(STATUS "CMAKE_INSTALL_LIBDIR = ${CMAKE_INSTALL_LIBDIR}")

# Optional compressor-based approximators (PDR decomposition strategies).
# The core library (multilevel decomposition, bitplane encoding, error control)
# only requires zstd; enable these to build the corresponding approximators
# and the tools that depend on them.
option(PRODM_WITH_SZ2 "Enable the SZ2 approximator" OFF)
option(PRODM_WITH_SZ3 "Enable the SZ3 approximator and SZ3-based tools" OFF)
option(PRODM_WITH_HPEZ "Enable the HPEZ (QoZ) approximator" OFF)
option(PRODM_WITH_MGARD "Enable the MGARD approximator" OFF)

# zstd (required): taken from the system installation
find_path(ZSTD_INCLUDE_DIR NAMES zstd.h
  HINTS /opt/homebrew/include /usr/local/include /usr/include
)
find_library(ZSTD_LIBRARY NAMES zstd libzstd
  HINTS /opt/homebrew/lib /usr/local/lib
)
if(NOT ZSTD_INCLUDE_DIR OR NOT ZSTD_LIBRARY)
  message(FATAL_ERROR
    "Could not find zstd: header=${ZSTD_INCLUDE_DIR}, lib=${ZSTD_LIBRARY}\n"
    "Please install libzstd (e.g., 'brew install zstd' or 'apt install libzstd-dev').")
endif()

add_library(${PROJECT_NAME} INTERFACE)

target_include_directories(${PROJECT_NAME} INTERFACE
  "${CMAKE_CURRENT_SOURCE_DIR}/include"
  "${ZSTD_INCLUDE_DIR}"
)

target_link_libraries(${PROJECT_NAME} INTERFACE
  "${ZSTD_LIBRARY}"
)

if(PRODM_WITH_SZ3)
  find_library(SZ3_LIBRARY NAMES SZ3c sz3c
    HINTS
      "${CMAKE_CURRENT_SOURCE_DIR}/external/SZ3/install"
      "${CMAKE_CURRENT_SOURCE_DIR}/external/SZ3/build"
    PATH_SUFFIXES lib lib64 src ${CMAKE_INSTALL_LIBDIR}
  )
  if(NOT SZ3_LIBRARY)
    message(FATAL_ERROR
      "PRODM_WITH_SZ3 is ON but the SZ3 library was not found under external/SZ3; "
      "build it first (see build_script.sh) or set SZ3_LIBRARY.")
  endif()
  set(SZ3_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/external/SZ3/install/include")
  target_include_directories(${PROJECT_NAME} INTERFACE "${SZ3_INCLUDE_DIR}")
  target_link_libraries(${PROJECT_NAME} INTERFACE "${SZ3_LIBRARY}")
  target_compile_definitions(${PROJECT_NAME} INTERFACE PRODM_HAVE_SZ3)
endif()

if(PRODM_WITH_SZ2)
  find_library(SZ2_LIBRARY NAMES SZ
    HINTS "${CMAKE_CURRENT_SOURCE_DIR}/external/SZ2/install"
    PATH_SUFFIXES lib lib64
  )
  if(NOT SZ2_LIBRARY)
    message(FATAL_ERROR
      "PRODM_WITH_SZ2 is ON but the SZ2 library was not found under external/SZ2; "
      "build it first (see build_script.sh) or set SZ2_LIBRARY.")
  endif()
  set(SZ2_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/external/SZ2/install/include")
  target_include_directories(${PROJECT_NAME} INTERFACE "${SZ2_INCLUDE_DIR}")
  target_link_libraries(${PROJECT_NAME} INTERFACE "${SZ2_LIBRARY}")
  target_compile_definitions(${PROJECT_NAME} INTERFACE PRODM_HAVE_SZ2)
endif()

if(PRODM_WITH_HPEZ)
  # QoZ (HPEZ) is used header-only
  set(QoZ_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/external/QoZ/install/include")
  if(NOT EXISTS "${QoZ_INCLUDE_DIR}/QoZ/api/sz.hpp")
    message(FATAL_ERROR
      "PRODM_WITH_HPEZ is ON but the QoZ headers were not found under external/QoZ; "
      "fetch them first (see build_script.sh) or set QoZ_INCLUDE_DIR.")
  endif()
  target_include_directories(${PROJECT_NAME} INTERFACE "${QoZ_INCLUDE_DIR}")
  target_compile_definitions(${PROJECT_NAME} INTERFACE PRODM_HAVE_HPEZ)
endif()

if(PRODM_WITH_MGARD)
  find_library(MGARD_LIBRARY NAMES mgard
    HINTS "${CMAKE_CURRENT_SOURCE_DIR}/external/MGARD/install-serial"
    PATH_SUFFIXES lib lib64 ${CMAKE_INSTALL_LIBDIR}
  )
  if(NOT MGARD_LIBRARY)
    message(FATAL_ERROR
      "PRODM_WITH_MGARD is ON but the MGARD library was not found under external/MGARD; "
      "build it first (see build_script.sh) or set MGARD_LIBRARY.")
  endif()
  set(MGARD_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/external/MGARD/install-serial/include")
  target_include_directories(${PROJECT_NAME} INTERFACE "${MGARD_INCLUDE_DIR}")
  target_link_libraries(${PROJECT_NAME} INTERFACE "${MGARD_LIBRARY}")
  target_compile_definitions(${PROJECT_NAME} INTERFACE PRODM_HAVE_MGARD)
endif()

install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/ DESTINATION include)

add_subdirectory(test)
add_subdirectory(prl_src)
