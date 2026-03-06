#ifndef UNIQUE_HEADER_FILE
#define UNIQUE_HEADER_FILE

#include <Eigen/Dense>
#include <map>
#include <vector>


void unique(const Eigen::MatrixXi& matrix,
            std::vector<int>& uniqueIndices,
            std::vector<int>& counts,
            std::vector<int>& inverse) {
  uniqueIndices.clear();
  counts.clear();
  inverse.clear();
  
  std::map<std::vector<int>, int> rowToIndex;
  std::vector<std::vector<int>> uniqueRows;
  
  for (int i = 0; i < matrix.rows(); ++i) {
    std::vector<int> row;
    row.reserve(matrix.cols());
    for (int j = 0; j < matrix.cols(); ++j) {
      row.push_back(matrix(i, j));
    }
    if (rowToIndex.find(row) == rowToIndex.end()) {
      rowToIndex[row] = uniqueIndices.size();
      uniqueIndices.push_back(i);
      uniqueRows.push_back(row);
    }
    inverse.push_back(rowToIndex[row]);
  }
  
  counts.resize(uniqueRows.size(), 0);
  for (int i = 0; i < matrix.rows(); ++i) {
    std::vector<int> row;
    row.reserve(matrix.cols());
    for (int j = 0; j < matrix.cols(); ++j) {
      row.push_back(matrix(i, j));  // Access each element directly
    }
    counts[rowToIndex[row]]++;
  }
}


#endif
