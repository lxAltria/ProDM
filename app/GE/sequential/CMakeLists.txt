# GE QoI tools: the GE synthesizer instantiates all approximators, so these
# need every optional compressor enabled.
# refactor_d64 and qoi_Vtot_d64 are built from test/ (shared with the Hurricane
# example); this directory adds the remaining GE-specific QoI retrievers.
if(PRODM_WITH_SZ2 AND PRODM_WITH_SZ3 AND PRODM_WITH_HPEZ AND PRODM_WITH_MGARD)
    function(add_my_executable target source)
        add_executable(${target} ${source})
        target_link_libraries(${target} ${PROJECT_NAME})
    endfunction()

    add_my_executable(qoi_T_d64 qoi_T_d64.cpp)
    add_my_executable(qoi_C_d64 qoi_C_d64.cpp)
    add_my_executable(qoi_Mach_d64 qoi_Mach_d64.cpp)
    add_my_executable(qoi_PT_d64 qoi_PT_d64.cpp)
    add_my_executable(qoi_mu_d64 qoi_mu_d64.cpp)
else()
    message(STATUS "Skipping the GE QoI tools (require PRODM_WITH_SZ2/SZ3/HPEZ/MGARD)")
endif()
