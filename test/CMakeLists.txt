function(add_my_executable target source)
    add_executable(${target} ${source})
    target_include_directories(${target} PRIVATE ${MGARDx_INCLUDES} ${SZ3_INCLUDES} ${ZSTD_INCLUDES})
    target_link_libraries(${target} ${PROJECT_NAME} ${SZ3_LIB} ${ZSTD_LIB})
endfunction()

add_my_executable(test_mdr_refactor test_mdr_refactor.cpp)
add_my_executable(test_mdr_reconstructor test_mdr_reconstructor.cpp)
add_my_executable(test_pdr_refactor test_pdr_refactor.cpp)
add_my_executable(test_pdr_reconstructor test_pdr_reconstructor.cpp)


