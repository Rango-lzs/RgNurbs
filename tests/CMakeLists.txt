cmake_minimum_required(VERSION 3.2)
project(rg_tests)

add_subdirectory(gtest_repo)

include_directories(${RG_HEADERS_DIR})
include_directories(gtest_repo/googletest/include)

#set(SOURCE_FILES main.cpp src/divider_tests.cpp)
file(GLOB SRC_FILES *.cpp *.h src/*.cpp src/*.h)

add_executable(rg_tests ${SRC_FILES})
target_link_libraries(rg_tests RgNurbs gtest)
install(TARGETS rg_tests DESTINATION ${RG_INSTALL_BIN_DIR})

