#include <cmath>
#include <iostream>
#include <vector>

#include "hungarian_algorithm.h"

int main() {
  // example 1
  std::vector<std::vector<double>> cost_matrix1 = {
      {8, 25, 50},
      {50, 35, 75},
      {22, 48, 150},
  };

  // example 2
  std::vector<std::vector<double>> cost_matrix2 = {
      {20, 15, 18, 20, 25}, {18, 20, 12, 14, 15}, {21, 23, 25, 27, 25},
      {17, 18, 21, 23, 20}, {18, 18, 16, 19, 20},
  };

  HungarianAlgorithm hungarian(cost_matrix1);
  HungarianSolution solution = hungarian.Solve();

  std::cout << "Solution : " << std::endl;
  for (auto pair : solution.pair) {
    std::cout << pair.first + 1 << ", " << pair.second + 1 << std::endl;
  }
  std::cout << std::endl;

  std::cout << "Cost : " << solution.cost << std::endl;

  return 0;
}
