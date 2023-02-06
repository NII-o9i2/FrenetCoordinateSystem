//
// Created by jojo-feng on 23-1-29.
//
#include "frenet_coordinate_system.h"
#include "frenet_coordinate_system_new.h"

#ifndef BENCHMARK_FRENET_COORDINATE_SYSTEM_FRENET_COORDINATE_MAIN_H_
#define BENCHMARK_FRENET_COORDINATE_SYSTEM_FRENET_COORDINATE_MAIN_H_
msquare::FrenetCoordinateSystem gen_frenet_coordinate(void);
void search_point(msquare::FrenetCoordinateSystem & frenet_coor);

msquare::FrenetCoordinateSystemNew gen_new_frenet_coordinate(void);
void search_point_new(msquare::FrenetCoordinateSystemNew & frenet_coor);
#endif //BENCHMARK_FRENET_COORDINATE_SYSTEM_FRENET_COORDINATE_MAIN_H_
