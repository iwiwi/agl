cmake_minimum_required(VERSION 2.8)
project(agl)

set(CMAKE_VERBOSE_MAKEFILE ON)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pthread -std=c++11")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -g")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g")

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

include_directories(. ../3rd_party)

file(GLOB TEST_SOURCE_FILES */*_test.cc */*/*_test.cc)
file(GLOB MAIN_SOURCE_FILES */*_main.cc */*/*_main.cc)
file(GLOB LIBRARY_SOURCE_FILES */*.h */*.cc)
list(REMOVE_ITEM LIBRARY_SOURCE_FILES ${TEST_SOURCE_FILES} ${MAIN_SOURCE_FILES})

foreach(MAIN_SOURCE_FILE IN LISTS MAIN_SOURCE_FILES)
    string(REGEX MATCH "([^/]+)_main.cc" NAME ${MAIN_SOURCE_FILE})
    set(RUNTIME_NAME ${CMAKE_MATCH_1})
    add_executable(${RUNTIME_NAME} ${LIBRARY_SOURCE_FILES} ${MAIN_SOURCE_FILE})
    target_link_libraries(${RUNTIME_NAME} gflags)
endforeach()

add_library(gtest STATIC
    ../3rd_party/gtest/gtest.h
    ../3rd_party/gtest/gtest_all.cc
    ../3rd_party/gtest/gtest_main.cc)

add_library(gflags STATIC
    ../3rd_party/gflags/config.h
    ../3rd_party/gflags/gflags_completions.cc
    ../3rd_party/gflags/gflags_completions.h
    ../3rd_party/gflags/gflags_declare.h
    ../3rd_party/gflags/gflags_gflags.h
    ../3rd_party/gflags/gflags_reporting.cc
    ../3rd_party/gflags/gflags.cc
    ../3rd_party/gflags/gflags.h
    ../3rd_party/gflags/mutex.h
    ../3rd_party/gflags/util.h)

add_executable(test_all EXCLUDE_FROM_ALL ${LIBRARY_SOURCE_FILES} ${TEST_SOURCE_FILES})
target_link_libraries(test_all gflags gtest)
