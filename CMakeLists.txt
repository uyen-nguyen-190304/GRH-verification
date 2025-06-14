cmake_minimum_required(VERSION 3.16)
project(grhverify LANGUAGES CXX)

# --------------------------------------------------------------------
# 1. Require pybind11 and Python3
# --------------------------------------------------------------------
find_package(pybind11 CONFIG REQUIRED)
find_package(Python3 COMPONENTS Interpreter Development REQUIRED)

# --------------------------------------------------------------------
# 2. Build the extension module
# --------------------------------------------------------------------
pybind11_add_module(
    grhverify                           # name of output *.so / *.pyd module
    MODULE                              # tells CMake it's a Python module
    src/rh_verification.cpp
    src/pybind.cpp
)

set_target_properties(grhverify PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/grhverify")
    
# --------------------------------------------------------------------
# 3. Use modern C++
# --------------------------------------------------------------------
target_compile_features(grhverify PRIVATE cxx_std_17)
if (MSVC)
    target_compile_options(grhverify PRIVATE /W4 /permissive-)
else()
    target_compile_options(grhverify PRIVATE -Wall -Wextra -Wpedantic)
endif()

# --------------------------------------------------------------------
# 4. Optional install support
# --------------------------------------------------------------------
install(TARGETS grhverify
        LIBRARY DESTINATION grhverify)  # used if someone runs `pip install -e .`
