set(TARGET_NAME RgNurbs)

project(${TARGET_NAME} LANGUAGES CXX)
set(SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})

if(MSVC)
    add_compile_options(
        $<$<CONFIG:Debug>:/MDd> #---|-- Statically link the runtime libraries
        $<$<CONFIG:Release>:/MD> #--|
    )
endif()

add_library(${TARGET_NAME} SHARED "")
target_compile_definitions(${TARGET_NAME} PRIVATE RG_NURBS_MODULE)
target_include_directories(${TARGET_NAME} PRIVATE ${SOURCE_DIR})

file(GLOB SRC_FILES *.cpp *.h)
source_group("" FILES ${SRC_FILES})
target_sources(${TARGET_NAME} PRIVATE ${SRC_FILES})

install(TARGETS ${TARGET_NAME} DESTINATION ${RG_INSTALL_BIN_DIR})