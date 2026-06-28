cmake_minimum_required(VERSION 3.23)
project(lenses LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(lenses INTERFACE)

target_include_directories(lenses
        INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_compile_features(lenses INTERFACE cxx_std_23)

add_executable(lenses_demo main.cpp)
target_link_libraries(lenses_demo PRIVATE lenses)

enable_testing()

add_executable(lenses_tests tests/trace_tests.cpp)
target_link_libraries(lenses_tests PRIVATE lenses)
add_test(NAME lenses_tests COMMAND lenses_tests)
