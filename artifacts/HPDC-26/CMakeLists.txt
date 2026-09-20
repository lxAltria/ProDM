# HPDC'26 (QProR): parallel weighted refactoring and QoI retrieval. parallel/para_refactor and para_Vtot (float, approximator
# selectable) are built; the older double-precision GE tools paraGE_* and para_refactor_d64 are kept as sources only.
# Needs MPI and all four approximators (the drivers go through the GE synthesizer recipes).
if(NOT PRODM_WITH_MPI)
    return()
endif()
find_package(MPI)
if(NOT MPI_FOUND)
    message(STATUS "Skipping the HPDC-26 parallel tools (MPI not found)")
    return()
endif()
if(PRODM_WITH_SZ2 AND PRODM_WITH_SZ3 AND PRODM_WITH_HPEZ AND PRODM_WITH_MGARD)
    add_executable(para_refactor parallel/para_refactor.cpp)
    target_link_libraries(para_refactor ${PROJECT_NAME} MPI::MPI_CXX)
    add_executable(para_Vtot parallel/para_Vtot.cpp)
    target_link_libraries(para_Vtot ${PROJECT_NAME} MPI::MPI_CXX)
else()
    message(STATUS "Skipping para_refactor and para_Vtot (require PRODM_WITH_SZ2/SZ3/HPEZ/MGARD)")
endif()
