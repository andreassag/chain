cmake_minimum_required(VERSION 3.20)
project(chain VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_compile_options(-Wall -Wextra -Wpedantic -O3)

include_directories(include)

add_executable(chain src/main.cpp)

enable_testing()
add_executable(chain_tests tests/test_chain.cpp)
add_test(NAME chain_tests COMMAND chain_tests)
