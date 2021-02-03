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

#include "Graph.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <locale>
#include <set>
#include <sstream>
#include <string>
#include <vector>
using std::ifstream;
using std::ofstream;
using std::set;
using std::string;
using std::stringstream;
using std::vector;
using std::endl;
using std::max;
using std::min;
using std::cerr;
using std::cout;

Graph::Graph(double modularity_resolution)
{
	m_size = 0;
	m_total_weight = 0.0;
	m_is_directed = false;
	m_number_of_communities = 0;
	m_modularity_resolution = modularity_resolution;
}


Graph::~Graph(void)
{
}

void Graph::FillMatrix(const vector<int>& sources, const vector<int>& destinations, const vector<double>& weights)
{
	if (!m_is_directed)
		m_total_weight *= 2;
	// we expect vertices to be numbered starting from 0
	m_size = 1 + max(*max_element(sources.begin(), sources.end()), *max_element(destinations.begin(), destinations.end()));
	m_matrix.assign(m_size, vector<double>(m_size, 0));
	for (size_t i = 0; i < sources.size(); ++i) {
		m_matrix[sources[i]][destinations[i]] += weights[i];
		if (!m_is_directed)
			m_matrix[destinations[i]][sources[i]] += weights[i];
	}
}

void Graph::FillModMatrix(const vector<int>& sources, const vector<int>& destinations, const vector<double>& weights)
{
	if (!m_is_directed)
		m_total_weight *= 2;
	m_modularity_matrix.assign(m_size, vector<double>(m_size, 0));
	vector<double> sumQ2(m_size, 0.0);
	vector<double> sumQ1(m_size, 0.0);
	for (size_t i = 0; i < sources.size(); ++i) {
		m_modularity_matrix[sources[i]][destinations[i]] += weights[i] / m_total_weight;
		if (!m_is_directed)
			m_modularity_matrix[destinations[i]][sources[i]] += weights[i] / m_total_weight;
		sumQ1[sources[i]] += weights[i] / m_total_weight;
		sumQ2[destinations[i]] += weights[i] / m_total_weight;
		if (!m_is_directed) {
			sumQ1[destinations[i]] += weights[i] / m_total_weight;
			sumQ2[sources[i]] += weights[i] / m_total_weight;
		}
	}
	for (int i = 0; i < m_size; ++i)
		for (int j = 0; j < m_size; ++j)
			m_modularity_matrix[i][j] -= m_modularity_resolution * sumQ1[i]*sumQ2[j];
	for (int i = 0; i < m_size; ++i)
		for (int j = 0; j < m_size; ++j)
			m_modularity_matrix[i][j] = m_modularity_matrix[j][i] = (m_modularity_matrix[i][j] + m_modularity_matrix[j][i]) / 2;
}

void Graph::ReadFromEdgelist(const string& file_name)
{
	ifstream file(file_name.c_str());
	if (!file.is_open()) {
        cerr << "File " << file_name << " can not be opened." << endl;
		return;
    }
	vector<int> sources, destinations;
	vector<double> weights;
	int min_vertex_number = 2e9;
	int max_vertex_number = 0;
	while (file.good()) {
		string line;
		getline(file, line);
		int src = -1, dst = -1;
		double weight = 1.0;
		stringstream str_stream(line);
		str_stream >> src >> dst;
		if (!str_stream.eof())
			str_stream >> weight;
		if (!str_stream.fail() && src != -1 && dst != -1)
		{
			min_vertex_number = min(min_vertex_number, min(src, dst));
			max_vertex_number = max(max_vertex_number, max(src, dst));
			sources.push_back(src);
			destinations.push_back(dst);
			weights.push_back(weight);
			m_total_weight += weight;
		}
	}
	file.close();
	m_is_directed = true;
	m_size = 1 + max_vertex_number - min_vertex_number;
	for (size_t i = 0; i < sources.size(); ++i) {
		sources[i] -= min_vertex_number;
		destinations[i] -= min_vertex_number;
	}
	FillModMatrix(sources, destinations, weights);
}

void Graph::ReadFromPajeck(const string& file_name)
{
	ifstream file(file_name.c_str());
	if (!file.is_open()) {
		cerr << "File " << file_name << " can not be opened." << endl;
		return;
    }
    std::locale locale;
	vector<int> sources, destinations;
	vector<double> weights;
	int min_vertex_number = 2e9;
	int max_vertex_number = 0;
	bool skip = true;
	while (file.good()) {
		string line;
		getline(file, line);
		stringstream str_stream(line);
        string trimmed_string;
        str_stream >> trimmed_string; // Strips carriage return on Windows
        transform(trimmed_string.begin(), trimmed_string.end(), trimmed_string.begin(), [&locale](char elem){return std::tolower(elem, locale);});
		if (trimmed_string == "*edges") {
			skip = false;
			m_is_directed = false;
		} else if (trimmed_string == "*arcs") {
			skip = false;
			m_is_directed = true;
		} else if (skip) {
			str_stream.str(line);
			int v = -1;
			str_stream >> v;
			if (!str_stream.fail()) {
				min_vertex_number = min(min_vertex_number, v);
				max_vertex_number = max(max_vertex_number, v);
			}
		} else if (!skip) {
			int src = -1, dst = -1;
			double weight = 1.0;
			str_stream.str(line);
			str_stream >> src >> dst;
			if (!str_stream.eof())
				str_stream >> weight;
			if (!str_stream.fail() && src != -1 && dst != -1) {
				min_vertex_number = min(min_vertex_number, min(src, dst));
				max_vertex_number = max(max_vertex_number, max(src, dst));
				sources.push_back(src);
				destinations.push_back(dst);
				weights.push_back(weight);
				m_total_weight += weight;
			}
		}
	}
	file.close();
	m_size = 1 + max_vertex_number - min_vertex_number;
	for (size_t i = 0; i < sources.size(); ++i) {
		sources[i] -= min_vertex_number;
		destinations[i] -= min_vertex_number;
	}
	FillModMatrix(sources, destinations, weights);
}

double Graph::EdgeWeight(int u, int v) const
{
	return m_matrix[u][v];
}

void Graph::CalcModMatrix()
{
	if (!m_modularity_matrix.empty())
		return;

	m_modularity_matrix.assign(m_size, vector<double>(m_size, 0.0));
	for (int i = 0; i < m_size; ++i)
		for (int j = 0; j < m_size; ++j)
			m_modularity_matrix[i][j] = EdgeWeight(i, j) / m_total_weight;
	
	vector<double> sumQ2(m_size, 0.0);
	vector<double> sumQ1(m_size, 0.0);
	for (int i = 0; i < m_size; ++i)
		for (int j = 0; j < m_size; ++j) {
			sumQ1[i] += m_modularity_matrix[i][j];
			sumQ2[j] += m_modularity_matrix[i][j];
		}
	for (int i = 0; i < m_size; ++i)
		for (int j = 0; j < m_size; ++j)
			m_modularity_matrix[i][j] -= m_modularity_resolution * sumQ1[i]*sumQ2[j];
	for (int i = 0; i < m_size; ++i)
		for (int j = 0; j < m_size; ++j)
			m_modularity_matrix[i][j] = m_modularity_matrix[j][i] = (m_modularity_matrix[i][j] + m_modularity_matrix[j][i]) / 2;
}

void Graph::Print() const
{
	cout << "Matrix:" << endl;
	for (int i = 0; i < m_size; ++i) {
		for (int j = 0; j < m_size; ++j) {
			cout << m_matrix[i][j] << '\t';
		}
		cout << endl;
	}
	cout << "Modularity matrix:" << endl;
	for (int i = 0; i < m_size; ++i) {
		for (int j = 0; j < m_size; ++j) {
			cout << m_modularity_matrix[i][j] << '\t';
		}
		cout << endl;
	}
}

void Graph::PrintCommunity(const string& file_name) const
{
	ofstream file(file_name.c_str());
	if (!file.is_open()) {
        cerr << "File " << file_name << " can not be opened." << endl;
		return;
    }
	for (int i = 0; i < m_size; ++i)
		file << m_communities[i] << endl;
	file.close();
}

void Graph::SetCommunities(const vector<int>& new_communities, int number)
{
	if (m_size != int(new_communities.size()))
		return;
	m_communities = new_communities;
	if (number == -1)
		m_number_of_communities = *max_element(m_communities.begin(), m_communities.end()) + 1;
	else
		m_number_of_communities = number;
}

double Graph::Modularity() const
{
	double modularity = 0;
	for (size_t i = 0; i < m_modularity_matrix.size(); ++i)
		for (size_t j = 0; j < m_modularity_matrix.size(); ++j)
			if (m_communities[i] == m_communities[j])
				modularity += m_modularity_matrix[i][j];
	return modularity;
}

void Graph::PerformSplit(int origin, int destination, const vector<int>& to_be_moved)
{
	if (destination > m_number_of_communities)
		destination = m_number_of_communities;
	if (destination == m_number_of_communities)
		++m_number_of_communities;
	for (int i = 0; i < m_size; ++i)
		if (m_communities[i] == origin && to_be_moved[i])
			m_communities[i] = destination;
}

bool Graph::IsCommunityEmpty(int community) const
{
	for (int i = 0; i < m_size; ++i)
		if (m_communities[i] == community)
			return false;
	return true;
}

bool Graph::DeleteCommunityIfEmpty(int community)
{
	if (IsCommunityEmpty(community)) {
		set<int> community_labels;
        for (int i = 0; i < m_size; ++i) {
			if (m_communities[i] > community)
				--m_communities[i];
			community_labels.insert(m_communities[i]);
		}
		m_number_of_communities = community_labels.size();
        return true;
	}
	return false;
}

vector<int> Graph::CommunityIndices(int community) const
{
	vector<int> res;
	for (int i = 0; i < m_size; ++i)
		if (m_communities[i] == community)
			res.push_back(i);
	return res;
}

vector< vector<double> > Graph::GetModularitySubmatrix(const vector<int>& indices) const
{
	vector< vector<double> > res(indices.size(), vector<double>(indices.size()));
	for (size_t i = 0; i < indices.size(); ++i)
		for (size_t j = 0; j < indices.size(); ++j)
			res[i][j] = m_modularity_matrix[indices[i]][indices[j]];
	return res;
}

vector<double> Graph::GetCorrectionVector(const vector<int>& orig_comm_ind, const vector<int>& dest_comm_ind) const
{
	vector<double> res(orig_comm_ind.size(), 0.0);
	for (size_t i = 0; i < orig_comm_ind.size(); ++i)
		for (size_t j = 0; j < dest_comm_ind.size(); ++j)
			res[i] += m_modularity_matrix[dest_comm_ind[j]][orig_comm_ind[i]];
	return res;
}
