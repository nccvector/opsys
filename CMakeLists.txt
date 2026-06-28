cmake_minimum_required(VERSION 3.23)
project(osys LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(osys INTERFACE)

target_include_directories(osys
        INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>)

target_sources(osys
        INTERFACE
        FILE_SET HEADERS
        BASE_DIRS include
        FILES
        include/osys/geometry.hpp
        include/osys/osys.hpp
        include/osys/medium.hpp
        include/osys/presets.hpp
        include/osys/optical_system.hpp
        include/osys/ray.hpp
        include/osys/sagitta.hpp)

target_compile_features(osys INTERFACE cxx_std_23)

add_executable(osys_demo main.cpp)
target_link_libraries(osys_demo PRIVATE osys)

find_package(raylib QUIET)
if (TARGET raylib)
        set(OSYS_RAYLIB_TARGET raylib)
else()
        find_package(PkgConfig QUIET)
        if (PkgConfig_FOUND)
                pkg_check_modules(RAYLIB QUIET IMPORTED_TARGET raylib)
                if (TARGET PkgConfig::RAYLIB)
                        set(OSYS_RAYLIB_TARGET PkgConfig::RAYLIB)
                endif()
        endif()
endif()

if (OSYS_RAYLIB_TARGET)
        add_executable(osys_editor osys_editor.cpp)
        target_link_libraries(osys_editor PRIVATE osys ${OSYS_RAYLIB_TARGET})
else()
        message(STATUS "raylib not found; osys_editor target disabled")
endif()

enable_testing()

add_executable(osys_tests tests/trace_tests.cpp)
target_link_libraries(osys_tests PRIVATE osys)
add_test(NAME osys_tests COMMAND osys_tests)
