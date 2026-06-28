cmake_minimum_required(VERSION 3.23)
project(lenses LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(lenses INTERFACE)

target_include_directories(lenses
        INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>)

target_sources(lenses
        INTERFACE
        FILE_SET HEADERS
        BASE_DIRS include
        FILES
        include/lenses/geometry.hpp
        include/lenses/lenses.hpp
        include/lenses/medium.hpp
        include/lenses/optical_system.hpp
        include/lenses/ray.hpp
        include/lenses/sagitta.hpp)

target_compile_features(lenses INTERFACE cxx_std_23)

add_executable(lenses_demo main.cpp)
target_link_libraries(lenses_demo PRIVATE lenses)

enable_testing()

add_executable(lenses_tests tests/trace_tests.cpp)
target_link_libraries(lenses_tests PRIVATE lenses)
add_test(NAME lenses_tests COMMAND lenses_tests)
