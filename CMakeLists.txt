cmake_minimum_required(VERSION 3.23)
project(opsys LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(opsys INTERFACE)

target_include_directories(opsys
        INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>)

target_sources(opsys
        INTERFACE
        FILE_SET HEADERS
        BASE_DIRS include
        FILES
        include/opsys/geometry.hpp
        include/opsys/opsys.hpp
        include/opsys/medium.hpp
        include/opsys/presets.hpp
        include/opsys/optical_system.hpp
        include/opsys/ray.hpp
        include/opsys/sagitta.hpp)

target_compile_features(opsys INTERFACE cxx_std_23)

add_executable(opsys_demo main.cpp)
target_link_libraries(opsys_demo PRIVATE opsys)

find_package(raylib QUIET)
if (TARGET raylib)
        set(OPSYS_RAYLIB_TARGET raylib)
else()
        find_package(PkgConfig QUIET)
        if (PkgConfig_FOUND)
                pkg_check_modules(RAYLIB QUIET IMPORTED_TARGET raylib)
                if (TARGET PkgConfig::RAYLIB)
                        set(OPSYS_RAYLIB_TARGET PkgConfig::RAYLIB)
                endif()
        endif()
endif()

if (OPSYS_RAYLIB_TARGET)
        add_executable(opsys_editor opsys_editor.cpp)
        target_link_libraries(opsys_editor PRIVATE opsys ${OPSYS_RAYLIB_TARGET})
else()
        message(STATUS "raylib not found; opsys_editor target disabled")
endif()

enable_testing()

add_executable(opsys_tests tests/trace_tests.cpp)
target_link_libraries(opsys_tests PRIVATE opsys)
add_test(NAME opsys_tests COMMAND opsys_tests)
