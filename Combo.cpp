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


#include "Combo.h"
#include "Graph.h"
#include "Matrix.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <numeric>
#include <optional>
#include <random>
#include <vector>
using namespace std;


double ModularityGain(const Matrix& Q, const vector<double>& correction_vector, const vector<int>& community)
{
	size_t n = community.size();
	double mod_gain = 0.0;
	for (size_t i = 0; i < n; ++i) {
		for (size_t j = 0; j < n; ++j)
			if (community[i] == community[j])
				mod_gain += Q[i][j];
			else
				mod_gain -= Q[i][j];
	}
	mod_gain *= 0.5;
	for (size_t i = 0; i < n; ++i) {
		if (community[i])
			mod_gain += correction_vector[i];
		else
			mod_gain -= correction_vector[i];
	}
	return mod_gain;
}

//perform a split improvement using a Kernighan-Lin-style iterative shifts series
double PerformKernighansShift(const Matrix& Q, const vector<double>& correction_vector,
	const vector<int>& communities_old, vector<int>& communities_new)
{
 	size_t n = Q.size();
	vector<double> gains(n, 0.0);
	for (size_t i = 0; i < n; ++i) {
		for (size_t j = 0; j < n; ++j)
			if (i != j) {
				if (communities_old[i] == communities_old[j])
					gains[i] -= Q[i][j];
				else
					gains[i] += Q[i][j];
			}
		if (communities_old[i])
			gains[i] -= correction_vector[i];
		else
			gains[i] += correction_vector[i];
		gains[i] *= 2;
	}
	vector<double> gains_got(n, 0.0);
	vector<size_t> gains_indexes(n, 0);
	communities_new = communities_old;
	for (size_t i = 0; i < n; ++i) {
		vector<double>::iterator it = max_element(gains.begin(), gains.end());
		gains_got[i] = *it;
		size_t gains_ind = size_t(it - gains.begin());
		gains_indexes[i] = gains_ind;
		if (i > 0)
			gains_got[i] = gains_got[i] + gains_got[i-1];
		for (size_t j = 0; j < n; ++j)
			if (communities_new[gains_ind] == communities_new[j])
				gains[j] += 4 * Q[gains_ind][j];
			else
				gains[j] -= 4 * Q[gains_ind][j];
		communities_new[gains_ind] = !communities_new[gains_ind];
		gains[gains_ind] = -INF;
	}
	vector<double>::iterator it = max_element(gains_got.begin(), gains_got.end());
	double mod_gain = *it;
	size_t steps_to_get_max_gain = size_t(it - gains_got.begin() + 1);
	if (mod_gain > 0) {
		communities_new = communities_old;
		for (size_t i = 0; i < steps_to_get_max_gain; ++i)
			communities_new[gains_indexes[i]] = !communities_new[gains_indexes[i]];
	} else {
		communities_new = communities_old;
		mod_gain = 0;
	}
	return mod_gain;
}

//try to split the subnetwork with respect to the correction vector
double ComboAlgorithm::Split(Matrix& Q,
	const vector<double>& correction_vector, vector<int>& to_be_moved)
{
	double mod_gain = 0.0;
	vector<double> sumQ = Sum(Q, 1);
	size_t n = Q.size();
	for (size_t i = 0; i < n; ++i)
		Q[i][i] += 2 * correction_vector[i] - sumQ[i]; //adjust the submatrix
	int tries;
	if (m_num_split_attempts > 0)
		tries = m_num_split_attempts;
	else
		tries = int(pow(abs(log(m_current_best_gain)), m_autoC2) / m_autoC1 + 3);
	for (int tryI = 1; tryI <= tries; ++tryI) {
		vector<int> communities(n); // 0 - stay in origin, 1 - move to destination
		//perform an initial simple split
		if (m_fixed_split_step > 0 && tryI <= 6 * m_fixed_split_step && tryI % m_fixed_split_step == 0) {
			//perform one of predefined split types
			int fixed_split_type = tryI / m_fixed_split_step;
			if (fixed_split_type == 1 || fixed_split_type == 2)
				communities.assign(n, 2 - fixed_split_type);
			else {
				vector<double> sum_pos = Sum(Q, 1, Positive);
				size_t node_ind;
				if (fixed_split_type == 3 || fixed_split_type == 4)
					node_ind = size_t(max_element(sum_pos.begin(), sum_pos.end()) - sum_pos.begin());
				else
					node_ind = size_t(min_element(sum_pos.begin(), sum_pos.end()) - sum_pos.begin());
				communities.assign(n, -1);
				int community = 1;
				communities[node_ind] = community;
				while (true) {
					optional<size_t> next_node_ind;
					double cur_min = 1e300;
					double cur_max = -1e300;
					for (size_t i = 0; i < n; ++i) {
						if (communities[i] == -1) {
							if ((fixed_split_type == 3 || fixed_split_type == 5) && Q[node_ind][i] < cur_min) {
								next_node_ind = i;
								cur_min = Q[node_ind][i];
							} else if ((fixed_split_type == 4 || fixed_split_type == 6) && Q[node_ind][i] > cur_max) {
								next_node_ind = i;
								cur_max = Q[node_ind][i];
							}
						}
					}
					if (!next_node_ind.has_value())
						break;
					node_ind = next_node_ind.value();
					community ^= 1;
					communities[node_ind] = community;
				}
			}
		} else {
			for (size_t i = 0; i < n; ++i)
				communities[i] = m_bernoulli_distribution(m_random_number_generator);
		}
		double mod_gain_total = ModularityGain(Q, correction_vector, communities);
		double mod_gain_from_shift = 1;
		while (mod_gain_from_shift > THRESHOLD) {
			vector<int> communities_shifted(n);
			mod_gain_from_shift = PerformKernighansShift(Q, correction_vector, communities, communities_shifted);
			if (mod_gain_from_shift > THRESHOLD) {
				mod_gain_total += mod_gain_from_shift;
				communities = communities_shifted;
				if (m_debug_verify) {
					double mod_gain_verify = ModularityGain(Q, correction_vector, communities);
					if (fabs(mod_gain_verify - mod_gain_total) > THRESHOLD)
						cerr << "ERROR" << endl;
				}
			}
		}
		if (mod_gain < mod_gain_total) {
			to_be_moved = communities;
			mod_gain = mod_gain_total;
		}
		if (mod_gain <= 1e-6)
			tries = int(tries / 2);
	}
	if (fabs(mod_gain) < THRESHOLD)
		to_be_moved.assign(n, 1);
	return mod_gain;
}

void ComboAlgorithm::ReCalc(Graph& graph, vector< vector<double> >& move_gains, vector< vector<bool> >& splits_communities, size_t origin, size_t destination)
{
	move_gains[origin][destination] = 0;
	if (origin != destination) {
		vector<size_t> orig_comm_ind = graph.CommunityIndices(origin);
		if (!orig_comm_ind.empty()) {
			vector<double> correction_vector = graph.GetCorrectionVector(orig_comm_ind, graph.CommunityIndices(destination));
			vector<int> to_be_moved(orig_comm_ind.size());
			Matrix Q = graph.GetModularitySubmatrix(orig_comm_ind);
			move_gains[origin][destination] = Split(Q, correction_vector, to_be_moved);
			for (size_t i = 0; i < to_be_moved.size(); ++i)
				splits_communities[destination][orig_comm_ind[i]] = to_be_moved[i];
		}
	}
}

double BestGain(const vector< vector<double> >& move_gains, size_t& origin, size_t& destination)
{
	double best_gain = -1;
	for (size_t i = 0; i < move_gains.size(); ++i)
		for (size_t j = 0; j < move_gains[i].size(); ++ j)
			if (best_gain < move_gains[i][j]) {
				best_gain = move_gains[i][j];
				origin = i;
				destination = j;
			}
	return best_gain;
}

bool DeleteCommunityIfEmpty(Graph& graph, vector< vector<double> >& move_gains, vector< vector<bool> >& splits_communities, size_t origin)
{
	if (graph.DeleteCommunityIfEmpty(origin)) {
		for (size_t i = origin; i+1 < move_gains.size(); ++i)
			move_gains[i] = move_gains[i+1];
		move_gains.back().assign(move_gains.back().size(), 0);
		for (size_t i = 0; i < move_gains.size(); ++i) {
			for (size_t j = origin; j+1 < move_gains[i].size(); ++j)
				move_gains[i][j] = move_gains[i][j+1];
			move_gains[i].back() = 0;
		}
		for (size_t i = origin; i+1 < splits_communities.size(); ++i)
			splits_communities[i] = splits_communities[i+1];
		splits_communities.back().assign(splits_communities.back().size(), false);
		return true;
	}
	return false;
}

void ComboAlgorithm::Run(Graph& graph, optional<size_t> max_communities, bool start_separate,
	optional<string> intermediate_result_file_name)
{
	if (!max_communities.has_value())
		max_communities = graph.Size();
	vector<size_t> initial_comm(graph.Size(), 0);
	if(start_separate)
		iota(initial_comm.begin(), initial_comm.end(), 0);
	graph.SetCommunities(initial_comm);
	double currentMod = graph.Modularity();
	if (m_output_info_level > 0) {
		cout << "0. " << graph.NumberOfCommunities() << " communities, "
			<< "initial modularity = " << currentMod << endl;
	}
	vector< vector<double> > move_gains(graph.NumberOfCommunities(),
		vector<double>(graph.NumberOfCommunities() + (graph.NumberOfCommunities() < max_communities), 0)); //results of splitting communities
	//vectors of boolean meaning that corresponding vertex should be moved to that destination
	vector< vector<bool> > splits_communities(graph.NumberOfCommunities() + (graph.NumberOfCommunities() < max_communities), vector<bool>(graph.Size(), false)); //best split vectors
	m_current_best_gain = 1;
	size_t origin = 0, destination = 0;
	for (origin = 0; origin < graph.NumberOfCommunities(); ++ origin)
		for (destination = 0; destination < graph.NumberOfCommunities() + (graph.NumberOfCommunities() < max_communities); ++destination)
			ReCalc(graph, move_gains, splits_communities, origin, destination);
	m_current_best_gain = BestGain(move_gains, origin, destination);
	int iteration = 0;
	while (m_current_best_gain > THRESHOLD) {
		++iteration;
		bool community_added = destination >= graph.NumberOfCommunities();
		graph.PerformSplit(origin, destination, splits_communities[destination]);
		bool origin_became_empty = DeleteCommunityIfEmpty(graph, move_gains, splits_communities, origin);
		if (origin_became_empty) {
			if (community_added)
				cerr << "WARNING: moving ALL nodes to EMPTY community should not occur." << endl;
			community_added = false;
			if (origin < destination)
				--destination;
		}
		if (m_output_info_level > 0) {
			cout << iteration << ". " << graph.NumberOfCommunities() << " communities, "
				<< "modularity = " << graph.Modularity() << ", last modularity gain = " << m_current_best_gain << endl;
		}
		if (intermediate_result_file_name.has_value() && intermediate_result_file_name.value() != "") {
			graph.PrintCommunity(intermediate_result_file_name.value());
		}
		if (m_debug_verify) {
			double oldMod = currentMod;
			currentMod = graph.Modularity();
			if (fabs(currentMod - oldMod - m_current_best_gain) > THRESHOLD)
				cerr << "ERROR: modularity does not match." << endl;
		}
		if (community_added) {
			if (destination + 1 < max_communities) {
				for (auto& row : move_gains)
					row.push_back(row[destination]);
				splits_communities.push_back(splits_communities[destination]);
			}
			if (destination >= move_gains.size())
				move_gains.push_back(vector<double>(move_gains.back().size(), 0));
		}
		for (size_t i = 0; i < graph.NumberOfCommunities() + (graph.NumberOfCommunities() < max_communities); ++i) {
			if (!origin_became_empty)
				ReCalc(graph, move_gains, splits_communities, origin, i);
			ReCalc(graph, move_gains, splits_communities, destination, i);
			if (!origin_became_empty && i != destination && i < graph.NumberOfCommunities())
				ReCalc(graph, move_gains, splits_communities, i, origin);
			if (i != origin && i < graph.NumberOfCommunities())
				ReCalc(graph, move_gains, splits_communities, i, destination);
		}
		m_current_best_gain = BestGain(move_gains, origin, destination);
	}
	if (m_output_info_level > 0) {
		cout << "Finished with " << graph.NumberOfCommunities() << " communities, "
			<< "achieved modularity = " << graph.Modularity() << endl;
	}
}

void ComboAlgorithm::SetNumberOfSplitAttempts(int split_tries)
{
	if (split_tries <= 0) {
		if (split_tries == -1) {
			m_autoC1 = 1.5*log(10);
			m_autoC2 = 1;
		} else if (split_tries == -2) {
			m_autoC1 = log(10);
			m_autoC2 = 1;
		} else {
            m_autoC1 = 2;
            m_autoC2 = 1.5;
		}  
	}
	m_num_split_attempts = split_tries;
}

ComboAlgorithm::ComboAlgorithm(optional<uint_fast32_t> random_seed, int num_split_attempts, int fixed_split_step, int output_info_level) :
	m_fixed_split_step(fixed_split_step),
	m_output_info_level(output_info_level),
	m_random_number_generator(random_seed.has_value() ? random_seed.value() :
		static_cast<uint_fast32_t>(std::chrono::duration_cast<std::chrono::microseconds>(
			std::chrono::steady_clock::now().time_since_epoch()).count())),
	m_bernoulli_distribution(0.5)
{
	SetNumberOfSplitAttempts(num_split_attempts);
}
