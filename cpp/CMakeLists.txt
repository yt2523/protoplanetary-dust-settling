cmake_minimum_required(VERSION 3.10)
project(protoplanetary_dust_settling_cpp LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# No fast-math; explicitly disable floating-point contraction for reproducibility.
# The flag is compiler-specific: GNU/Clang accept -ffp-contract=off, MSVC uses /fp:strict.
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  add_compile_options(-ffp-contract=off)
elseif(MSVC)
  add_compile_options(/fp:strict)
endif()

add_library(dust_settling_lib STATIC src/dust_settling.cpp)
target_include_directories(dust_settling_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)

add_executable(dust_settling src/main.cpp)
target_link_libraries(dust_settling PRIVATE dust_settling_lib)
