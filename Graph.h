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

#include <string>
#include <vector>

class Graph
{
public:
	explicit Graph(double modularity_resolution = 1);
	virtual ~Graph();

	void ReadFromEdgelist(const std::string& file_name);
	void ReadFromPajeck(const std::string& file_name);
	void CalcModMatrix();

	int Size() const {return m_size;}
	int NumberOfCommunities() const {return m_number_of_communities;};
	double EdgeWeight(int u, int v) const;
	bool IsCommunityEmpty(int community) const;

	double Modularity() const;
	std::vector< std::vector<double> > GetModularitySubmatrix(const std::vector<int>& indices) const;
	std::vector<double> GetCorrectionVector(const std::vector<int>& orig_comm_ind, const std::vector<int>& dest_comm_ind) const;
	
	void SetCommunities(const std::vector<int>& new_communities, int number = -1);
	std::vector<int> Communities() const {return m_communities;};
	std::vector<int> CommunityIndices(int comm) const;

	void PerformSplit(int origin, int destination, const std::vector<int>& split_communities);
	bool DeleteCommunityIfEmpty(int community);
	void Print() const;
	void PrintCommunity(const std::string& file_name) const;

private:
	void FillMatrix(const std::vector<int>& sources, const std::vector<int>& destinations, const std::vector<double>& weights);
	void FillModMatrix(const std::vector<int>& sources, const std::vector<int>& destinations, const std::vector<double>& weights);

private:
	int m_size;
	double m_total_weight;
	int m_number_of_communities;
	bool m_is_directed;
	// Modularity Resolution Parameter
	// as per Newman 2016 (https://journals.aps.org/pre/abstract/10.1103/PhysRevE.94.052315)
	double m_modularity_resolution;
	std::vector<std::vector<double> > m_matrix;
	std::vector<std::vector<double> > m_modularity_matrix;
	std::vector<int> m_communities;
};

#endif //GRAPH_H
