#include <cmath>
#include <iostream>
#include <vector>

#include "hungarian_algorithm.h"

struct Point {
  Point(double x, double y) : x(x), y(y) {};
  double x, y;
};

int main() {
  std::vector<Point> io_pad = {
      Point(15000, 45000), Point(20000, 45000), Point(25000, 45000),
      Point(30000, 45000), Point(35000, 45000), Point(45000, 35000),
      Point(45000, 30000), Point(45000, 25000), Point(45000, 20000),
      Point(45000, 15000), Point(35000, 5000),  Point(30000, 5000),
      Point(25000, 5000),  Point(20000, 5000),  Point(15000, 5000),
      Point(5000, 15000),  Point(5000, 20000),  Point(5000, 25000),
      Point(5000, 30000),  Point(5000, 35000)};

  std::vector<Point> rdl = {
      Point(15000, 32500), Point(20000, 32500), Point(25000, 32500),
      Point(30000, 32500), Point(35000, 32500), Point(15000, 27500),
      Point(20000, 27500), Point(25000, 27500), Point(30000, 27500),
      Point(35000, 27500), Point(15000, 22500), Point(20000, 22500),
      Point(25000, 22500), Point(30000, 22500), Point(35000, 22500),
      Point(15000, 17500), Point(20000, 17500), Point(25000, 17500),
      Point(30000, 17500), Point(35000, 17500)};

  std::vector<std::vector<double>> cost_matrix;

  for (int i = 0; i < 20; i++) {
    cost_matrix.push_back(std::vector<double>());
    for (int j = 0; j < 20; j++) {
      double value =
          std::abs(io_pad[i].x - rdl[j].x) + std::abs(io_pad[i].y - rdl[j].y);
      cost_matrix[i].push_back(value);
    }
  }

  // std::vector<std::vector<double>> cost_matrix = {
  //     {20, 15, 18, 20, 25}, {18, 20, 12, 14, 15}, {21, 23, 25, 27, 25},
  //     {17, 18, 21, 23, 20}, {18, 18, 16, 19, 20},
  // };

  HungarianAlgorithm hungarian(cost_matrix);
  hungarian.Solve();

  return 0;
}
