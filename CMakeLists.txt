cmake_minimum_required(VERSION 3.10)
project(cell_growth_simulation)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add include directories
include_directories(${PROJECT_SOURCE_DIR}/include)

# Add source files
set(SOURCES
    src/main.cpp
    src/cell_simulation.cpp
    src/particle.cpp
)

# Create executable
add_executable(cell_simulation ${SOURCES})

# Add compiler warnings
if(MSVC)
    target_compile_options(cell_simulation PRIVATE /W4)
else()
    target_compile_options(cell_simulation PRIVATE -Wall -Wextra -Wpedantic)
endif() 