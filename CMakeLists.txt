cmake_minimum_required(VERSION 3.2)

project(RgNurbs)

set(CMAKE_CXX_STANDARD 17)

if(MSVC)
    add_compile_options(
        $<$<CONFIG:Debug>:/MDd> #---|-- Statically link the runtime libraries
        $<$<CONFIG:Release>:/MD> #--|
    )
endif()

message(STATUS "This info is for debug: " ${CMAKE_BINARY_DIR})

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)
set_property(GLOBAL PROPERTY PREDEFINED_TARGETS_FOLDER "CMakePredefinedTargets")


set(CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR})

set(RG_INSTALL_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)
set(RG_INSTALL_BIN_DIR ${PROJECT_SOURCE_DIR}/bin)
set(RG_INSTALL_LIB_DIR ${PROJECT_SOURCE_DIR}/lib)

set(RG_HEADERS_DIR ${PROJECT_SOURCE_DIR}/src)

include_directories(${RG_INSTALL_INCLUDE_DIR})
include_directories(${RG_HEADERS_DIR})

add_subdirectory(src)
add_subdirectory(tests)

set_target_properties(
    gtest gtest_main gmock gmock_main
    PROPERTIES FOLDER "gtest"
)
