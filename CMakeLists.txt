cmake_minimum_required(VERSION 3.16)
project(grhverify LANGUAGES CXX)

# --------------------------------------------------------------------
# 1.  Locate or download pybind11
# --------------------------------------------------------------------
include(FetchContent)

# First try a system / pip install
find_package(pybind11 CONFIG QUIET)

# If not found, fetch from GitHub (one-time download, then cached)
if(NOT pybind11_FOUND)
    message(STATUS "pybind11 not found – fetching with FetchContent")
    FetchContent_Declare(
        pybind11
        GIT_REPOSITORY https://github.com/pybind/pybind11.git
        GIT_TAG        v2.11.1         # ← pin a release
    )
    FetchContent_MakeAvailable(pybind11)
endif()

# --------------------------------------------------------------------
# 2.  Locate Python  (for include & linker flags)
# --------------------------------------------------------------------
find_package(Python3 COMPONENTS Interpreter Development REQUIRED)

# --------------------------------------------------------------------
# 3.  Build the pybind11 extension
# --------------------------------------------------------------------
pybind11_add_module(
    grhverify                # output name  grhverify.cpython-*.so
    MODULE
    src/rh_verification.cpp
    src/pybind.cpp
)

# place the shared library directly into grhverify/ so that
# `import grhverify` works without PYTHONPATH gymnastics
set_target_properties(grhverify PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/grhverify")

# --------------------------------------------------------------------
# 4.  Compile flags
# --------------------------------------------------------------------
target_compile_features(grhverify PRIVATE cxx_std_17)

if(MSVC)
    target_compile_options(grhverify PRIVATE /W4 /permissive-)
else()
    target_compile_options(grhverify PRIVATE -Wall -Wextra -Wpedantic)
    # optimise if we are in Release (default)
    if(CMAKE_BUILD_TYPE STREQUAL "" OR CMAKE_BUILD_TYPE STREQUAL "Release")
        target_compile_options(grhverify PRIVATE -O3 -march=native)
    endif()
endif()

# --------------------------------------------------------------------
# 5.  (Optional) install target – useful for  `pip install -e .`
# --------------------------------------------------------------------
install(TARGETS grhverify
        LIBRARY DESTINATION grhverify)
