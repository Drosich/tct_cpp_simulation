cmake_minimum_required(VERSION 3.15)
project(tct_sim LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


# --- Find ROOT ---
# This assumes you have sourced ROOT’s environment:
#   source /path/to/root/bin/thisroot.sh
find_package(ROOT REQUIRED COMPONENTS RIO Net Hist Graf Graf3d Gpad)
find_package(nlohmann_json REQUIRED)
include(${ROOT_USE_FILE})
# Add include directory
include_directories(include
                    ${ROOT_INCLUDE_DIRS}
                    /usr/include)

# Collect all source files from src/
file(GLOB_RECURSE SRC_FILES src/*.cc)

# Add the executable (main.cpp + all source files)
add_executable(${PROJECT_NAME}
    main.cc
    ${SRC_FILES}
)

# --- Link against ROOT libraries ---
target_link_libraries(${PROJECT_NAME} PUBLIC ${ROOT_LIBRARIES} nlohmann_json::nlohmann_json)
target_include_directories(${PROJECT_NAME} PUBLIC ${ROOT_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/include /usr/include)


# Optional: extra warnings (for GCC/Clang)
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(${PROJECT_NAME} PRIVATE -Wall -Wextra -pedantic)
endif()

