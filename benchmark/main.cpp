#include <iostream>
#include "frenet_coordinate_system/frenet_coordinate_main.h"

#include <benchmark/benchmark.h>

static void BM_StringCreation(benchmark::State& state) {
  for (auto _ : state)
    std::string empty_string;
}
// Register the function as a benchmark
BENCHMARK(BM_StringCreation);

// Define another benchmark
static void BM_StringCopy(benchmark::State& state) {
  std::string x = "hello";
  for (auto _ : state)
    std::string copy(x);
}
BENCHMARK(BM_StringCopy);

// frenet update
static void BM_Frenet_Update(benchmark::State& state){
//  auto frenet_coor = gen_frenet_coordinate();
  for (auto _ : state){
    auto frenet_coor = gen_frenet_coordinate();
//    search_point(frenet_coor);
  }
}
BENCHMARK(BM_Frenet_Update);

// frenet search points
static void BM_Frenet_Search(benchmark::State& state){
  auto frenet_coor = gen_frenet_coordinate();
  for (auto _ : state){
    search_point(frenet_coor);
  }
}
BENCHMARK(BM_Frenet_Search);


// new frenet update
static void BM_New_Frenet_Update(benchmark::State& state){
// auto frenet_coor = gen_new_frenet_coordinate();
 for (auto _ : state){
  auto frenet_coor = gen_new_frenet_coordinate();
 }
}
BENCHMARK(BM_New_Frenet_Update);

// new frenet search
static void BM_New_Frenet_Search(benchmark::State& state){
 auto frenet_coor = gen_new_frenet_coordinate();
  for (auto _ : state){
    search_point_new(frenet_coor);
  }
}
BENCHMARK(BM_New_Frenet_Search);


// hypot test
static void BM_hypot_double(benchmark::State& state){
  double double_x = 5.0;
  double double_y = 7.0;
  for (auto _ : state){
    for (int i = 0; i<10000;i++) {
      auto tmp = std::hypotf(double_x*i, double_y*i);
    }
  }
}
BENCHMARK(BM_hypot_double);


// hypotf test
static void BM_hypotf(benchmark::State& state){
  double double_x = 5.0;
  double double_y = 7.0;
  for (auto _ : state){
    for (int i = 0; i<10000;i++) {
      auto tmp = std::hypotf(double_x*i, double_y*i);
    }
  }
}
BENCHMARK(BM_hypotf);


// hypot test
static void BM_hypot(benchmark::State& state){
  double double_x = 5.0;
  double double_y = 7.0;
  for (auto _ : state){
    for (int i = 0; i<10000;i++){
      auto tmp = std::hypot(double_x*i,double_y*i);
    }
  }
}
BENCHMARK(BM_hypot);

BENCHMARK_MAIN();


