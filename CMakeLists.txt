cmake_minimum_required(VERSION 3.13)
project(CQOI VERSION 0.0.0.1 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_FLAGS "-g -O3")

include(GNUInstallDirs)
message(STATUS "CMAKE_INSTALL_LIBDIR = ${CMAKE_INSTALL_LIBDIR}")

find_library(SZ2_LIBRARY NAMES SZ
  HINTS "${CMAKE_CURRENT_SOURCE_DIR}/external/SZ2/install"
  PATH_SUFFIXES lib lib64
)

find_package(ADIOS2 REQUIRED CONFIG
  HINTS "${CMAKE_CURRENT_SOURCE_DIR}/external/ADIOS2/adios2-install"
  PATH_SUFFIXES ${CMAKE_INSTALL_LIBDIR}/cmake/adios2
)

find_path(ZSTD_INCLUDE_DIR NAMES zstd.h
  HINTS
    "${CMAKE_CURRENT_SOURCE_DIR}/external/SZ2/install/include"
    "${CMAKE_CURRENT_SOURCE_DIR}/external/SZ3/tools/zstd/lib"
    "${CMAKE_CURRENT_SOURCE_DIR}/external/SZ3/tools/zstd/lib/common"
    /usr/include /usr/local/include
)

find_library(ZSTD_LIBRARY NAMES zstd libzstd
  HINTS "${CMAKE_CURRENT_SOURCE_DIR}/external/SZ2/install"
  PATH_SUFFIXES lib ${CMAKE_INSTALL_LIBDIR}
)

if(NOT ZSTD_INCLUDE_DIR OR NOT ZSTD_LIBRARY)
  message(FATAL_ERROR
    "Could not find zstd: header=${ZSTD_INCLUDE_DIR}, lib=${ZSTD_LIBRARY}\n"
    "Please ensure libzstd is installed or built in SZ2/SZ3.")
endif()

find_library(SZ3_LIBRARY NAMES SZ3c sz3c
  HINTS
    "${CMAKE_CURRENT_SOURCE_DIR}/external/SZ3/install"
    "${CMAKE_CURRENT_SOURCE_DIR}/external/SZ3/build"
  PATH_SUFFIXES lib lib64 src ${CMAKE_INSTALL_LIBDIR}
)

find_library(MGARD_LIBRARY NAMES mgard
  HINTS "${CMAKE_CURRENT_SOURCE_DIR}/external/MGARD/install-serial"
  PATH_SUFFIXES lib lib64 ${CMAKE_INSTALL_LIBDIR}
)

if(NOT MGARD_LIBRARY)
  message(FATAL_ERROR "Could not find MGARD library")
endif()

set(SZ3_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/external/SZ3/install/include")

if(NOT SZ3_LIBRARY)
  message(FATAL_ERROR
    "Could not find SZ3 library.\n"
    "Tried paths:\n"
    "  ${CMAKE_CURRENT_SOURCE_DIR}/external/SZ3/install/lib/\n"
    "  ${CMAKE_CURRENT_SOURCE_DIR}/external/SZ3/install/lib64/\n"
    "  ${CMAKE_CURRENT_SOURCE_DIR}/external/SZ3/build/src/\n"
    "  ${CMAKE_CURRENT_SOURCE_DIR}/external/SZ3/build/lib/\n"
  )
endif()

# find_package(QoZ REQUIRED CONFIG
#   HINTS "${CMAKE_CURRENT_SOURCE_DIR}/external/QoZ/install" 
# )

set(QoZ_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/external/QoZ/install/include")

set(SZ2_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/external/SZ2/install/include")

set(MGARD_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/external/MGARD/install-serial/include")

message(STATUS "QoZ_INCLUDE_DIR = ${QoZ_INCLUDE_DIR}")

set(MGARDx_INCLUDES "${CMAKE_CURRENT_SOURCE_DIR}/external/MGARDx/install/include")

add_library(${PROJECT_NAME} INTERFACE)

target_include_directories(${PROJECT_NAME} INTERFACE
  "${CMAKE_CURRENT_SOURCE_DIR}/include"
  "${SZ3_INCLUDE_DIR}"
  "${ZSTD_INCLUDE_DIR}"
  "${QoZ_INCLUDE_DIR}"
  "${SZ2_INCLUDE_DIR}"
  "${MGARD_INCLUDE_DIR}"
)

target_link_libraries(${PROJECT_NAME} INTERFACE
  "${ZSTD_LIBRARY}"
  "${SZ3_LIBRARY}"
  "${SZ2_LIBRARY}"
  "${MGARD_LIBRARY}"
)

install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/ DESTINATION include)

add_subdirectory(test)
add_subdirectory(prl_src)