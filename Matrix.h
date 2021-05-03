/*
 Copyright 2021
 Alexander Belyi <alexander.belyi@gmail.com>,
 Stanislav Sobolevsky <stanly@mit.edu>
 
 This file is part of Combo algorithm.
 
 Combo is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 Combo is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with Combo.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef MATRIX_H
#define MATRIX_H

#include <cstdlib>
#include <iostream>
#include <limits>
#include <set>
#include <utility>
#include <vector>

const double INF = std::numeric_limits<double>::max();
const double THRESHOLD = 1e-6;
const double EPS = 2.0 * std::numeric_limits<double>::epsilon();

typedef std::vector<std::vector<double>> Matrix;
typedef std::vector<std::vector<int>> MatrixInt;
typedef std::vector<std::vector<bool>> MatrixBool;

template<typename T> bool Positive(T x) {return x > 0.0;}
template<typename T> bool Negative(T x) {return x < 0.0;}
template<typename T> bool NotNegative(T x) {return x >= 0.0;}
template<typename T> bool NotPositive(T x) {return x <= 0.0;}
template<typename T> std::vector<T> Sum(const std::vector<std::vector<T>>& matrix, int dim, bool (*Pred)(T) = NULL)
{
    if (dim == 1) {
        std::vector<T> res(matrix.size(), 0);
        for (size_t i = 0; i < matrix.size(); ++i)
            for (size_t j = 0; j < matrix[i].size(); ++j)
                if (Pred == NULL)
                    res[i] = res[i] + matrix[i][j];
                else if (Pred(matrix[i][j]))
                    res[i] = res[i] + matrix[i][j];
        return res;
    } else {
        std::vector<T> res(matrix[0].size(), 0);
        for (size_t i = 0; i < matrix.size(); ++i)
            for (size_t j = 0; j < matrix[i].size(); ++j)
                if (Pred == NULL)
                    res[j] = res[j] + matrix[i][j];
                else if (Pred(matrix[i][j]))
                    res[j] = res[j] + matrix[i][j];
        return res;
    }
}

template<typename T> T Sum(const std::vector<T>& vec, bool (*Pred)(T) = NULL)
{
    size_t n = vec.size();
    T res = 0;
    for (size_t i = 0; i < n; ++i)
        if (Pred == NULL || Pred(vec[i]))
            res += vec[i];
    return res;
}

template<typename T> bool TestAll(const std::vector<T>& vec, bool (*Pred)(T))
{
    for (T element: vec)
        if (!Pred(element))
            return false;
    return true;
}

template<typename T> std::vector<T> diag(const std::vector<std::vector<T>>& matrix)
{
    size_t n = matrix.size();
    std::vector<T> res(n, 0.0);
    for (size_t i = 0; i < n; ++i)
        res[i] = matrix[i][i];
    return res;
}

template<typename T> std::vector<std::vector<T>> Submatrix(const std::vector<std::vector<T>>& matrix,
                                                           const std::set<size_t>& indices)
{
    size_t n = indices.size();
    std::vector<std::vector<T>> res(n, std::vector<T>(n));
    size_t i = 0, j = 0;
    for (size_t ind1 : indices) {
        j = 0;
        for (size_t ind2 : indices) {
            res[i][j] = matrix[ind1][ind2];
            ++j;
        }
        ++i;
    }
    return res;
}

template<typename T> std::vector<std::vector<T>> Submatrix(const std::vector<std::vector<T>>& matrix,
                                                           const std::vector<size_t>& indices)
{
    size_t n = indices.size();
    std::vector<std::vector<T>> res(n, std::vector<T>(n));
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < n; ++j)
            res[i][j] = matrix[indices[i]][indices[j]];
    return res;
}

template<typename T> std::vector<std::vector<T>> MergeTwoNodes(const std::vector<std::vector<T>>& matrix,
                                                               size_t node1, size_t node2)
{
    if (node1 > node2)
        std::swap(node1, node2);
    if (node1 == node2)
        return matrix;
    size_t n = matrix.size();
    std::vector<std::vector<T>> merged_matrix(n - 1, std::vector<T>(n - 1));
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            size_t merged_i = i, merged_j = j;
            if (i == node2)
                merged_i = node1;
            if (j == node2)
                merged_j = node1;
            if (merged_i > node2)
                --merged_i;
            if (merged_j > node2)
                --merged_j;
            merged_matrix[merged_i][merged_j] += matrix[i][j];
        }
    }
    return merged_matrix;
}

template<typename T>
T GetEdgeNumber(T i, T j, T n)
{
    /* edge numbering example for 5x5 adj matrix:
        0 1 2 3 4
      -----------
     0| - 0 1 2 3
     1| 0 - 4 5 6
     2| 1 4 - 7 8
     3| 2 5 7 - 9
     4| 3 6 8 9 -
     */
    if (j == i)
        std::cerr << "ERROR: we do not consider loops!" << std::endl;
    if (j > i)
        std::swap(i, j);
    return n * j + i - (j + 2) * (j + 1) / 2;
}

#endif //MATRIX_H
