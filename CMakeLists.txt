cmake_minimum_required(VERSION 3.14)

project(math)

add_library(math INTERFACE)
target_include_directories(math PUBLIC INTERFACE include)

add_subdirectory(test)
