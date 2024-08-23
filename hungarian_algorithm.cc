#include "hungarian_algorithm.h"

#include <iostream>

HungarianAlgorithm::HungarianAlgorithm(
    std::vector<std::vector<double>> cost_matrix)
    : source_matrix(cost_matrix), cost_matrix(cost_matrix), row(0), col(0) {
  row = cost_matrix.size();
  col = cost_matrix[0].size();
}

HungarianSolution HungarianAlgorithm::Solve() {
  int iteration = 0;

  while (true) {
#ifdef DEBUG
    std::cout << "Iteration " << iteration << std::endl;
    std::cout << std::endl;
#endif

    SubtractMinRowValue();

#ifdef DEBUG
    std::cout << "SubtractMinRowValue : " << std::endl;
    for (int i = 0; i < row; i++) {
      for (int j = 0; j < col; j++) {
        std::cout << cost_matrix[i][j] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
#endif

    SubtractMinColValue();

#ifdef DEBUG
    std::cout << "SubtractMinColValue : " << std::endl;
    for (int i = 0; i < row; i++) {
      for (int j = 0; j < col; j++) {
        std::cout << cost_matrix[i][j] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
#endif

    std::shared_ptr<bool> cover_row(new bool[row], [](bool* p) { delete[] p; });
    std::shared_ptr<bool> cover_col(new bool[row], [](bool* p) { delete[] p; });
    std::shared_ptr<short> element_type(new short[row * col],
                                        [](short* p) { delete[] p; });

    FindMinimumLines(cover_row, cover_col, element_type);

#ifdef DEBUG
    std::cout << "FindMinimumLines : " << std::endl;
    std::cout << std::endl;
    std::cout << "Cover Row : " << std::endl;
    for (int i = 0; i < row; i++) {
      std::cout << i << " : " << cover_row.get()[i] << std::endl;
    }
    std::cout << std::endl;

    std::cout << "Cover Col :" << std::endl;
    for (int i = 0; i < row; i++) {
      std::cout << i << " : " << cover_col.get()[i] << std::endl;
    }
    std::cout << std::endl;

    std::cout << "Element Type :" << std::endl;
    for (int i = 0; i < row; i++) {
      for (int j = 0; j < col; j++) {
        std::cout << element_type.get()[i * col + j] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
#endif

    int line_count = TotalLine(cover_row, cover_col);

    if (line_count != row) {
      double min = FindMinimumNumber(cover_row, cover_col);
      AdjustWeight(cover_row, cover_col, min);
    } else {
#ifdef DEBUG
      std::cout << "Find the optimal solution" << std::endl << std::endl;
#endif
      break;
    }

#ifdef DEBUG
    std::cout << "AdjustWeight : " << std::endl;
    for (int i = 0; i < row; i++) {
      for (int j = 0; j < col; j++) {
        std::cout << cost_matrix[i][j] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
#endif

    iteration++;
  }

  std::list<std::pair<int, int>> pair;
  std::shared_ptr<bool> choice_table(new bool[row],
                                     [](bool* p) { delete[] p; });
  std::fill(choice_table.get(), choice_table.get() + row, false);
  bool find = false;
  int cost = 0;

  FindOptimalSolution(pair, choice_table, cost, 0, find);

  return HungarianSolution(pair, cost);
}

void HungarianAlgorithm::SubtractMinRowValue() {
  for (int i = 0; i < row; i++) {
    double min = std::numeric_limits<double>::max();

    for (int j = 0; j < col; j++) {
      if (cost_matrix[i][j] < min) {
        min = cost_matrix[i][j];
      }
    }

    for (int j = 0; j < col; j++) {
      cost_matrix[i][j] -= min;
    }
  }
}

void HungarianAlgorithm::SubtractMinColValue() {
  std::unique_ptr<double[]> min_array(new double[col]);
  std::fill(min_array.get(), min_array.get(),
            std::numeric_limits<double>::max());

  for (int i = 0; i < row; i++) {
    for (int j = 0; j < col; j++) {
      if (cost_matrix[i][j] < min_array[j]) {
        min_array[j] = cost_matrix[i][j];
      }
    }
  }

  for (int i = 0; i < row; i++) {
    for (int j = 0; j < col; j++) {
      cost_matrix[i][j] -= min_array[j];
    }
  }
}

// References : https://en.wikipedia.org/wiki/Hungarian_algorithm
void HungarianAlgorithm::FindMinimumLines(std::shared_ptr<bool> cover_row,
                                          std::shared_ptr<bool> cover_col,
                                          std::shared_ptr<short> element_type) {
  std::fill(cover_row.get(), cover_row.get() + row, false);
  std::fill(cover_col.get(), cover_col.get() + row, false);
  std::fill(element_type.get(), element_type.get() + row * col, 0);

  for (int i = 0; i < row; i++) {
    for (int j = 0; j < col; j++) {
      if (cost_matrix[i][j] == 0) {
        if (!cover_col.get()[j]) {
          // star 0
          element_type.get()[i * col + j] = 1;
          cover_col.get()[j] = true;
          break;
        }
      }
    }
  }

cover_all_star:
  std::fill(cover_row.get(), cover_row.get() + row, false);
  std::fill(cover_col.get(), cover_col.get() + row, false);
  RemovePrime(element_type);

  short* element_type_pointer = element_type.get();

  for (int i = 0; i < row; i++) {
    for (int j = 0; j < col; j++) {
      if (element_type_pointer[i * col + j] == 1) {
        cover_col.get()[j] = true;
      }
    }
  }

find_non_covered_zero:
  for (int i = 0; i < row; i++) {
    for (int j = 0; j < col; j++) {
      if (cost_matrix[i][j] == 0 && element_type.get()[i * col + j] == 0 &&
          !cover_row.get()[i] && !cover_col.get()[j]) {
        // prime 0
        int index = FindStarInRow(element_type, i);
        element_type.get()[i * col + j] = 2;

        if (index != -1) {
          cover_row.get()[i] = true;
          cover_col.get()[index] = false;
          goto find_non_covered_zero;
        } else {
          std::shared_ptr<short> path_martix(new short[row * col],
                                             [](short* p) { delete[] p; });
          std::copy(element_type.get(), element_type.get() + row * col,
                    path_martix.get());
          int row_index = i;
          int col_index = j;
          element_type.get()[i * col + j] = 1;

          while (true) {
            row_index = FindStarInCol(path_martix, col_index);
            if (row_index != -1) {
              element_type.get()[row_index * col + col_index] = 0;
              col_index = FindPrimeInRow(path_martix, row_index);
              element_type.get()[row_index * col + col_index] = 1;
            } else {
              break;
            }
          }
          goto cover_all_star;
        }
      }
    }
  }
}

int HungarianAlgorithm::FindStarInRow(std::shared_ptr<short> element_type,
                                      int index) {
  for (int i = 0; i < col; i++) {
    if (element_type.get()[index * col + i] == 1) {
      return i;
    }
  }

  return -1;
}

int HungarianAlgorithm::FindStarInCol(std::shared_ptr<short> path_martix,
                                      int index) {
  for (int i = 0; i < col; i++) {
    if (path_martix.get()[i * col + index] == 1) {
      return i;
    }
  }
  return -1;
}

int HungarianAlgorithm::FindPrimeInRow(std::shared_ptr<short> path_martix,
                                       int index) {
  for (int i = 0; i < col; i++) {
    if (path_martix.get()[index * col + i] == 2) {
      return i;
    }
  }

  return -1;
}

void HungarianAlgorithm::RemovePrime(std::shared_ptr<short> element_type) {
  for (int i = 0; i < row; i++) {
    for (int j = 0; j < col; j++) {
      if (element_type.get()[i * col + j] == 2) {
        element_type.get()[i * col + j] = 0;
      }
    }
  }
}

int HungarianAlgorithm::TotalLine(std::shared_ptr<bool> cover_row,
                                  std::shared_ptr<bool> cover_col) {
  int count = 0;

  for (int i = 0; i < row; i++) {
    if (cover_row.get()[i]) {
      count++;
    }
    if (cover_col.get()[i]) {
      count++;
    }
  }

  return count;
}

double HungarianAlgorithm::FindMinimumNumber(std::shared_ptr<bool> cover_row,
                                             std::shared_ptr<bool> cover_col) {
  double min = std::numeric_limits<double>::max();

  for (int i = 0; i < row; i++) {
    for (int j = 0; j < col; j++) {
      if (cost_matrix[i][j] < min && !cover_row.get()[i] &&
          !cover_col.get()[j]) {
        min = cost_matrix[i][j];
      }
    }
  }

  return min;
}

void HungarianAlgorithm::AdjustWeight(std::shared_ptr<bool> cover_row,
                                      std::shared_ptr<bool> cover_col,
                                      double value) {
  for (int i = 0; i < row; i++) {
    for (int j = 0; j < col; j++) {
      cost_matrix[i][j] -= value;

      if (cover_row.get()[i] || cover_col.get()[j]) {
        cost_matrix[i][j] += value;
      }

      if (cover_row.get()[i] && cover_col.get()[j]) {
        cost_matrix[i][j] += value;
      }
    }
  }
}

void HungarianAlgorithm::FindOptimalSolution(
    std::list<std::pair<int, int>>& pair, std::shared_ptr<bool> choice_table,
    int& cost, int deep, bool& find) {
  for (int i = 0; i < col; i++) {
    if (cost_matrix[deep][i] == 0 && !choice_table.get()[i]) {
      pair.emplace_back(std::pair<int, int>(deep, i));
      choice_table.get()[i] = true;
      cost += source_matrix[deep][i];

      if (deep == row - 1) {
        find = true;
        return;
      }

      FindOptimalSolution(pair, choice_table, cost, deep + 1, find);

      if (find) {
        return;
      } else {
        pair.pop_back();
        choice_table.get()[i] = false;
        cost -= source_matrix[deep][i];
      }
    }
  }
}
