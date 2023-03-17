# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

file(MAKE_DIRECTORY
  "/home/jojo-feng/Downloads/benchmark/cmake-build-debug/_deps/osqp-src"
  "/home/jojo-feng/Downloads/benchmark/cmake-build-debug/_deps/osqp-build"
  "/home/jojo-feng/Downloads/benchmark/cmake-build-debug/_deps/osqp-subbuild/osqp-populate-prefix"
  "/home/jojo-feng/Downloads/benchmark/cmake-build-debug/_deps/osqp-subbuild/osqp-populate-prefix/tmp"
  "/home/jojo-feng/Downloads/benchmark/cmake-build-debug/_deps/osqp-subbuild/osqp-populate-prefix/src/osqp-populate-stamp"
  "/home/jojo-feng/Downloads/benchmark/cmake-build-debug/_deps/osqp-subbuild/osqp-populate-prefix/src"
  "/home/jojo-feng/Downloads/benchmark/cmake-build-debug/_deps/osqp-subbuild/osqp-populate-prefix/src/osqp-populate-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/home/jojo-feng/Downloads/benchmark/cmake-build-debug/_deps/osqp-subbuild/osqp-populate-prefix/src/osqp-populate-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/home/jojo-feng/Downloads/benchmark/cmake-build-debug/_deps/osqp-subbuild/osqp-populate-prefix/src/osqp-populate-stamp${cfgdir}") # cfgdir has leading slash
endif()
