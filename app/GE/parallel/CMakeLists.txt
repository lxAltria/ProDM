find_package(MPI REQUIRED)

function(add_my_executable_prl target source)
    add_executable(${target} ${source})
    target_include_directories(${target} PRIVATE ${MDR_INCLUDES} ${MGARDx_INCLUDES} ${SZ3_INCLUDES} ${ZSTD_INCLUDES} ${MPI_INCLUDE_PATH})
    target_link_libraries(${target} ${PROJECT_NAME} adios2::adios2 ${SZ3_LIB} ${ZSTD_LIB} ${MPI_C_LIBRARIES})
endfunction()

add_my_executable_prl(refactorGE refactorGE.cpp)
add_my_executable_prl(refactorGE_SZ3 refactorGE_SZ3.cpp)
add_my_executable_prl(refactorGE_SZ3_delta refactorGE_SZ3_delta.cpp)
add_my_executable_prl(testGE_VTOT testGE_VTOT.cpp)
add_my_executable_prl(testGE_VTOT_SZ3 testGE_VTOT_SZ3.cpp)
add_my_executable_prl(testGE_VTOT_SZ3_delta testGE_VTOT_SZ3_delta.cpp)