set(TARGET_NAME RgNurbs)

project(${TARGET_NAME} LANGUAGES CXX)
set(SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})

add_library(${TARGET_NAME} SHARED "")
target_compile_definitions(${TARGET_NAME} PRIVATE RG_NURBS_MODULE)

target_include_directories(${TARGET_NAME} PRIVATE ${SOURCE_DIR})

file(GLOB rootfiles *.cpp *.h)

source_group("" FILES ${rootfiles})
target_sources(${TARGET_NAME} PRIVATE ${rootfiles})
SUBDIRLIST(SUBDIRS ${SOURCE_DIR})
foreach(subdir ${SUBDIRS})
    file(GLOB subdirFiles ${subdir}/*.cpp ${subdir}/*.h)
    string(REPLACE "/" "\\" subdir ${subdir})
    source_group(${subdir} FILES ${subdirFiles})
    target_sources(${TARGET_NAME} PRIVATE ${subdirFiles})
endforeach()