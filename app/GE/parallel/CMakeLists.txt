# Parallel GE tools: require MPI, ADIOS2, and the SZ3 baselines.
# Skipped gracefully when any prerequisite is missing 
if(NOT PRODM_WITH_SZ3)
    message(STATUS "Skipping the parallel GE tools (require PRODM_WITH_SZ3)")
    return()
endif()
find_package(MPI)
if(NOT MPI_FOUND)
    message(STATUS "Skipping the parallel GE tools (MPI not found)")
    return()
endif()
find_package(ADIOS2 CONFIG QUIET
  HINTS "${CMAKE_SOURCE_DIR}/external/ADIOS2/adios2-install"
)
if(NOT ADIOS2_FOUND)
    message(STATUS "Skipping the parallel GE tools (ADIOS2 not found; see build_script.sh)")
    return()
endif()

function(add_my_executable_prl target source)
    add_executable(${target} ${source})
    target_link_libraries(${target} ${PROJECT_NAME} MPI::MPI_CXX adios2::adios2)
endfunction()

add_my_executable_prl(refactorGE refactorGE.cpp)
add_my_executable_prl(refactorGE_SZ3 refactorGE_SZ3.cpp)
add_my_executable_prl(refactorGE_SZ3_delta refactorGE_SZ3_delta.cpp)
add_my_executable_prl(testGE_VTOT testGE_VTOT.cpp)
add_my_executable_prl(testGE_VTOT_SZ3 testGE_VTOT_SZ3.cpp)
add_my_executable_prl(testGE_VTOT_SZ3_delta testGE_VTOT_SZ3_delta.cpp)
