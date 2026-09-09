function(add_my_executable target source)
    add_executable(${target} ${source})
    target_link_libraries(${target} ${PROJECT_NAME})
endfunction()

# core tools: multilevel decomposition + bitplane encoding, no optional dependencies
add_my_executable(test_mdr_refactor test_mdr_refactor.cpp)
add_my_executable(test_mdr_reconstructor test_mdr_reconstructor.cpp)
add_my_executable(test_proaicd_refactor test_proaicd_refactor.cpp)
add_my_executable(test_proaicd_reconstructor test_proaicd_reconstructor.cpp)

# PDR tools: the Dummy approximator always works; other approximators require
# the corresponding PRODM_WITH_* options
add_my_executable(test_pdr_refactor test_pdr_refactor.cpp)
add_my_executable(test_pdr_reconstructor test_pdr_reconstructor.cpp)
add_my_executable(test_pdr_delta_refactor test_pdr_delta_refactor.cpp)
add_my_executable(test_pdr_delta_reconstructor test_pdr_delta_reconstructor.cpp)

# QoI tools with HPEZ-based refactoring only
if(PRODM_WITH_HPEZ)
    add_my_executable(test_qoi_refactor test_qoi_refactor.cpp)
    add_my_executable(test_qoi_reconstructor test_qoi_reconstructor.cpp)
else()
    message(STATUS "Skipping test_qoi_refactor and test_qoi_reconstructor (require PRODM_WITH_HPEZ)")
endif()
