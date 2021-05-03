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

#ifndef GRAPH_H
#define GRAPH_H

#include "Matrix.h"

#include <string>
#include <tuple>
#include <vector>

class Graph;

Graph ReadGraphFromFile(const std::string& file_name, double mod_resolution = 1, bool treat_as_modularity = false);

class Graph
{
public:
	explicit Graph(bool is_directed = false, double modularity_resolution = 1);
	Graph(size_t size, const std::vector<int>& sources, const std::vector<int>& destinations, const std::vector<double>& weights,
		bool is_directed, double modularity_resolution = 1, bool treat_as_modularity = false);
	Graph(size_t size, const std::vector<std::tuple<int, int, double>>& edges, bool is_directed,
		double modularity_resolution = 1, bool treat_as_modularity = false);
	explicit Graph(const std::vector<std::vector<double>>& matrix,
		double modularity_resolution = 1, bool treat_as_modularity = false);
	explicit Graph(std::vector<std::vector<double>>&& matrix,
		double modularity_resolution = 1, bool treat_as_modularity = false);
	Graph(const Graph& graph);
	Graph(Graph&& graph) noexcept;
	Graph& operator=(Graph graph);

	size_t Size() const {return m_modularity_matrix.size();}
	bool IsDirected() const {return m_is_directed;}
	double ModularityResolution() const {return m_modularity_resolution;}
	size_t NumberOfCommunities() const {return m_number_of_communities;}

	double Modularity() const;
	Matrix GetModularityMatrix() const {return m_modularity_matrix;}
	Matrix GetModularitySubmatrix(const std::vector<size_t>& indices) const;
	std::vector<double> GetCorrectionVector(const std::vector<size_t>& orig_comm_ind, const std::vector<size_t>& dest_comm_ind) const;
	
	void SetCommunities(const std::vector<size_t>& new_communities, size_t number = 0);
	std::vector<size_t> Communities() const {return m_communities;}
	std::vector<size_t> CommunityIndices(size_t comm) const;
	bool IsCommunityEmpty(size_t community) const;

	void PerformSplit(size_t origin, size_t destination, const std::vector<bool>& split_communities);
	bool DeleteCommunityIfEmpty(size_t community);
	void Print() const;
	void PrintCommunity(const std::string& file_name = "") const;

	friend void swap(Graph& left, Graph& right);

protected:
	void CalcModMatrix(size_t size, const std::vector<int>& sources, const std::vector<int>& destinations, const std::vector<double>& weights);
	void FillModMatrix(size_t size, const std::vector<int>& sources, const std::vector<int>& destinations, const std::vector<double>& weights);
	void CalcModMatrix(size_t size, const std::vector<std::tuple<int, int, double>>& edges);
	void FillModMatrix(size_t size, const std::vector<std::tuple<int, int, double>>& edges);
	void CalcModMatrix(const Matrix& matrix);
	void FillModMatrix(const Matrix& matrix);
	void FillModMatrix(Matrix&& matrix);

private:
	size_t m_number_of_communities;
	bool m_is_directed;
	// Modularity Resolution Parameter
	// as per Newman 2016 (https://journals.aps.org/pre/abstract/10.1103/PhysRevE.94.052315)
	double m_modularity_resolution;
	Matrix m_modularity_matrix;
	std::vector<size_t> m_communities;
};

#endif //GRAPH_H
