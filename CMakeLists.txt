cmake_minimum_required(VERSION 3.14)

project(math)

set(MATH_SRC 
    src/cylindrical.cpp
    src/euler.cpp
    src/matrix2x2.cpp
    src/matrix3x3.cpp
    src/matrix4x4.cpp
    src/polar.cpp
    src/quaternion.cpp
    src/spherical.cpp
    src/vector2.cpp
    src/vector3.cpp
    src/vector4.cpp
)

add_library(math STATIC ${MATH_SRC})
target_include_directories(math PUBLIC include)

add_subdirectory(test)
