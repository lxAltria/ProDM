function(add_my_executable target source)
    add_executable(${target} ${source})
    target_link_libraries(${target} ${PROJECT_NAME})
endfunction()

# prodm_refactor / prodm_retrieve: the command line entry points of all pipelines
# (mdr, proaicd, pdr, pdr-delta; the approximators follow the PRODM_WITH_* options).
# example/test_script.sh exercises them.
add_my_executable(prodm_refactor prodm_refactor.cpp)
add_my_executable(prodm_retrieve prodm_retrieve.cpp)
