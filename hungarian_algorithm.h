#pragma once

#include <algorithm>
#include <limits>
#include <list>
#include <memory>
#include <utility>
#include <vector>

struct HungarianSolution {
  HungarianSolution(std::list<std::pair<int, int>> pair, int cost)
      : pair(pair), cost(cost) {};

  std::list<std::pair<int, int>> pair;
  int cost;
};

class HungarianAlgorithm {
 public:
  HungarianAlgorithm(std::vector<std::vector<double>> cost_matrix);

  HungarianSolution Solve();

 private:
  void SubtractMinRowValue();
  void SubtractMinColValue();
  void FindMinimumLines(std::shared_ptr<bool> cover_row,
                        std::shared_ptr<bool> cover_col,
                        std::shared_ptr<short> element_type);
  int FindStarInRow(std::shared_ptr<short> element_type, int index);
  int FindStarInCol(std::shared_ptr<short> path_martix, int index);
  int FindPrimeInRow(std::shared_ptr<short> path_martix, int index);
  void RemovePrime(std::shared_ptr<short> element_type);
  int TotalLine(std::shared_ptr<bool> cover_row,
                std::shared_ptr<bool> cover_col);
  double FindMinimumNumber(std::shared_ptr<bool> cover_row,
                           std::shared_ptr<bool> cover_col);
  void AdjustWeight(std::shared_ptr<bool> cover_row,
                    std::shared_ptr<bool> cover_col, double value);
  void FindOptimalSolution(std::list<std::pair<int, int>>& solution,
                           std::shared_ptr<bool> choice_table, int& cost,
                           int deep, bool& find);

  std::vector<std::vector<double>> source_matrix;
  std::vector<std::vector<double>> cost_matrix;
  int row;
  int col;
};
