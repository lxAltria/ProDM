find_package(MPI)
if(NOT MPI_FOUND)
    message(STATUS "MPI not found; skipping the parallel tools in prl_src (install MPI or set PRODM_WITH_MPI=OFF to silence)")
    return()
endif()

function(add_my_executable target source)
    add_executable(${target} ${source})
    target_link_libraries(${target} ${PROJECT_NAME} MPI::MPI_CXX)
endfunction()

# core parallel tools
add_my_executable(para_two_modes_refactor para_two_modes_refactor.cpp)
add_my_executable(para_two_modes_reconstructor para_two_modes_reconstructor.cpp)
add_my_executable(para_mdr_refactor para_mdr_refactor.cpp)
add_my_executable(para_mdr_reconstructor para_mdr_reconstructor.cpp)

# SZ3-delta baselines
if(PRODM_WITH_SZ3)
    add_my_executable(para_PSZ3-delta_refactor para_PSZ3-delta_refactor.cpp)
    add_my_executable(para_PSZ3-delta_reconstructor para_PSZ3-delta_reconstructor.cpp)
else()
    message(STATUS "Skipping para_PSZ3-delta tools (require PRODM_WITH_SZ3)")
endif()

# QoI tools: the GE synthesizer instantiates all approximators
if(PRODM_WITH_SZ2 AND PRODM_WITH_SZ3 AND PRODM_WITH_HPEZ AND PRODM_WITH_MGARD)
    add_my_executable(para_refactor para_refactor.cpp)
    add_my_executable(para_Vtot para_Vtot.cpp)
else()
    message(STATUS "Skipping para_refactor and para_Vtot (require PRODM_WITH_SZ2/SZ3/HPEZ/MGARD)")
endif()

# add_my_executable(refactorROD refactorROD.cpp)
# add_my_executable(test_PT test_PT.cpp)
# add_my_executable(para_refactor_d64 para_refactor_d64.cpp)
# add_my_executable(paraGE_Vtot paraGE_Vtot.cpp)
# add_my_executable(paraGE_T paraGE_T.cpp)
# add_my_executable(paraGE_C paraGE_C.cpp)
# add_my_executable(paraGE_Mach paraGE_Mach.cpp)
# add_my_executable(paraGE_PT paraGE_PT.cpp)
# add_my_executable(paraGE_mu paraGE_mu.cpp)
