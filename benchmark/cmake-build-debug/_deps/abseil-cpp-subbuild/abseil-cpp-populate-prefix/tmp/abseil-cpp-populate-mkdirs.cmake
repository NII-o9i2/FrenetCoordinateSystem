# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

file(MAKE_DIRECTORY
  "/home/jojo-feng/Downloads/benchmark/cmake-build-debug/_deps/abseil-cpp-src"
  "/home/jojo-feng/Downloads/benchmark/cmake-build-debug/_deps/abseil-cpp-build"
  "/home/jojo-feng/Downloads/benchmark/cmake-build-debug/_deps/abseil-cpp-subbuild/abseil-cpp-populate-prefix"
  "/home/jojo-feng/Downloads/benchmark/cmake-build-debug/_deps/abseil-cpp-subbuild/abseil-cpp-populate-prefix/tmp"
  "/home/jojo-feng/Downloads/benchmark/cmake-build-debug/_deps/abseil-cpp-subbuild/abseil-cpp-populate-prefix/src/abseil-cpp-populate-stamp"
  "/home/jojo-feng/Downloads/benchmark/cmake-build-debug/_deps/abseil-cpp-subbuild/abseil-cpp-populate-prefix/src"
  "/home/jojo-feng/Downloads/benchmark/cmake-build-debug/_deps/abseil-cpp-subbuild/abseil-cpp-populate-prefix/src/abseil-cpp-populate-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/home/jojo-feng/Downloads/benchmark/cmake-build-debug/_deps/abseil-cpp-subbuild/abseil-cpp-populate-prefix/src/abseil-cpp-populate-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/home/jojo-feng/Downloads/benchmark/cmake-build-debug/_deps/abseil-cpp-subbuild/abseil-cpp-populate-prefix/src/abseil-cpp-populate-stamp${cfgdir}") # cfgdir has leading slash
endif()
