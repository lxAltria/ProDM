function(add_my_executable target source)
    add_executable(${target} ${source})
    target_link_libraries(${target} ${PROJECT_NAME})
endfunction()

# core tools: multilevel decomposition + bitplane encoding, no optional dependencies
add_my_executable(test_mdr_refactor test_mdr_refactor.cpp)
add_my_executable(test_mdr_reconstructor test_mdr_reconstructor.cpp)
add_my_executable(two_modes_refactor two_modes_refactor.cpp)
add_my_executable(two_modes_reconstructor two_modes_reconstructor.cpp)

# PDR tools: the Dummy approximator always works; other approximators require
# the corresponding PRODM_WITH_* options
add_my_executable(test_pdr_refactor test_pdr_refactor.cpp)
add_my_executable(test_pdr_reconstructor test_pdr_reconstructor.cpp)
add_my_executable(test_pdr_delta_refactor test_pdr_delta_refactor.cpp)
add_my_executable(test_pdr_delta_reconstructor test_pdr_delta_reconstructor.cpp)

# QoI tools: the GE synthesizer instantiates all approximators, so these need
# every optional dependency enabled
if(PRODM_WITH_SZ2 AND PRODM_WITH_SZ3 AND PRODM_WITH_HPEZ AND PRODM_WITH_MGARD)
    add_my_executable(refactor_d64 refactor_d64.cpp)
    add_my_executable(qoi_Vtot_d64 qoi_Vtot_d64.cpp)
else()
    message(STATUS "Skipping refactor_d64 and qoi_Vtot_d64 (require PRODM_WITH_SZ2/SZ3/HPEZ/MGARD)")
endif()
