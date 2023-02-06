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



BENCHMARK_MAIN();


