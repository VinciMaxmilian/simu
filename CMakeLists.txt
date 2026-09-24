cmake_minimum_required(VERSION 3.16)
project(simu CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native")

find_package(OpenMP REQUIRED)

add_executable(simu
  src/main.cpp
  src/core/research.cpp
  src/core/emergence.cpp
  src/core/biochemistry.cpp
  src/labs/toy.cpp
  src/labs/chemistry.cpp
  src/labs/materials.cpp
  src/labs/nuclear.cpp
)
target_include_directories(simu PRIVATE src)
target_compile_options(simu PRIVATE -Wall -Wextra)
target_link_libraries(simu PRIVATE OpenMP::OpenMP_CXX)
