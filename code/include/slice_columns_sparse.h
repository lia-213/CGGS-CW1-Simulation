#ifndef SLICE_COLUMNS_SPARSE_HEADER_FILE
#define SLICE_COLUMNS_SPARSE_HEADER_FILE

#include <Eigen/Dense>
#include <Eigen/Sparse>


Eigen::SparseMatrix<double> slice_columns_sparse(const Eigen::SparseMatrix<double>& mat,
                                                 const Eigen::VectorXi& colIndices)
{
    int numRows = mat.rows();
    
    Eigen::SparseMatrix<double> result(numRows, colIndices.size());
    std::vector<Eigen::Triplet<double>> triplets;
    
    for (size_t i = 0; i < colIndices.size(); ++i) {
        int col = colIndices(i);
        for (Eigen::SparseMatrix<double>::InnerIterator it(mat, col); it; ++it) {
            triplets.emplace_back(it.row(), i, it.value());
        }
    }
    
    result.setFromTriplets(triplets.begin(), triplets.end());
    return result;
}


#endif
