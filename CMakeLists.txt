############################################################
# Copyright (c) 2026
# Agnieszka Sorensen
############################################################

cmake_minimum_required(VERSION 3.10)

project(SMASH_ROOT_Analysis VERSION 1.0 LANGUAGES CXX)



# Set the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)


########################################
# External libraries and software
########################################
# Find ROOT
find_package(ROOT REQUIRED)
# Load ROOT macros for CMake
include(${ROOT_USE_FILE})


########################################
# Copy shell script for copying files
########################################
configure_file(
    ${CMAKE_SOURCE_DIR}/input/copy_src_and_input_files_to_build.sh
    ${CMAKE_BINARY_DIR}/copy_src_and_input_files_to_build.sh
    COPYONLY
)



########################################
# Compile SMASH_ROOT_Analysis
########################################
# Add source files
set(SRC_FILES
./src/main.cc

./src/basic_functions.cc
./src/config.cc
./src/flow.cc
./src/multiplicity.cc
./src/plots.cc
./src/read_Particles.cc
./src/SMASH_config_info.cc
./src/yields.cc
)

# Add header files
set(HEADER_FILES
./src/basic_functions.h
./src/centrality.h
./src/config.h
./src/constants.h
./src/flow.h
./src/fourvector.h
./src/multiplicity.h
./src/plots.h
./src/range.h
./src/read_Particles.h
./src/SMASH_config_info.h
./src/yields.h
)

# Create the executable
add_executable(SMASH_ROOT_Analysis ${SRC_FILES})

# Include your headers
target_include_directories(SMASH_ROOT_Analysis
    PRIVATE
    ${PROJECT_SOURCE_DIR}/src
)

# Link libraries
target_link_libraries(SMASH_ROOT_Analysis
    PRIVATE
    ${ROOT_LIBRARIES}
)

# Set high-performance compiler flags
target_compile_options(SMASH_ROOT_Analysis PRIVATE -O3 -march=native)



########################################
# Copy example config to build directory
########################################
add_custom_command(
    TARGET SMASH_ROOT_Analysis POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
            ${CMAKE_SOURCE_DIR}/input/example_analysis_config.txt
            $<TARGET_FILE_DIR:SMASH_ROOT_Analysis>/analysis_config.txt
)


