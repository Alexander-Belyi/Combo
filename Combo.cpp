/*                                                                            
    Copyright 2021
    Alexander Belyi <alexander.belyi@gmail.com>,
    Stanislav Sobolevsky <stanly@mit.edu>                                               
                                                                            
    This is the main file of Combo algorithm.

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

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <random>
using namespace std;

#define THRESHOLD 1e-6


vector<double> Sum(const vector< vector<double> >& matrix)
{
	int n = matrix.size();
	vector<double> res(n, 0.0);
	for(int i = 0; i < n; ++i)
		for(int j = 0; j < n; ++j)
			res[i] += matrix[i][j];
	return res;
}

template<typename T> bool Positive(T x) {return x > 0.0;}
template<typename T> bool Negative(T x) {return x < 0.0;}
template<typename T> bool NotNegative(T x) {return x >= 0.0;}
template<typename T> bool NotPositive(T x) {return x <= 0.0;}
vector<double> SumPos(const vector< vector<double> >& matrix, bool (*Pred)(double) = NULL)
{
	int n = matrix.size();
	vector<double> res(n, 0.0);
	for(int i = 0; i < n; ++i)
		for(int j = 0; j < n; ++j)
			if(Pred && Pred(matrix[i][j]))
				res[i] += matrix[i][j];
	return res;
}

template<typename T>
bool TestAll(const vector<T>& vec, bool (*Pred)(T))
{
	int n = vec.size();
	for(int i = 0; i < n; ++i)
		if(!Pred(vec[i]))
			return false;
	return true;
}

double ModGain(const vector< vector<double> >& Q, const vector<double>& correctionVector, const vector<int>& community)
{
	int n = community.size();
	double mod_gain = 0.0;
	for(int i = 0; i < n; ++i)
	{
		for(int j = 0; j < n; ++j)
			if(community[i] == community[j])
				mod_gain += Q[i][j];
			else
				mod_gain -= Q[i][j];
	}
	mod_gain *= 0.5;
	for(int i = 0; i < n; ++i)
	{
		if(community[i])
			mod_gain += correctionVector[i];
		else
			mod_gain -= correctionVector[i];
	}
	return mod_gain;
}

double PerformKernighansShift(const vector< vector<double> >& Q, const vector<double>& correctionVector, const vector<int>& communitiesOld, vector<int>& communitiesNew) //perform a split improvement using a Karnigan-Lin-style iterative shifts series
{
 	int n = Q.size();
	vector<double> gains(n, 0.0);
	for(int i = 0; i < n; ++i)
	{
		for(int j = 0; j < n; ++j)
			if(i != j) {
				if(communitiesOld[i] == communitiesOld[j])
					gains[i] -= Q[i][j];
				else
					gains[i] += Q[i][j];
			}
		if(communitiesOld[i])
			gains[i] -= correctionVector[i];
		else
			gains[i] += correctionVector[i];
		gains[i] *= 2;
	}
	vector<double> gains_got(n, 0.0);
	vector<int> gains_indexes(n, 0);
	communitiesNew = communitiesOld;
	for(int i = 0; i < n; ++i)
	{
		vector<double>::iterator it = max_element(gains.begin(), gains.end());
		gains_got[i] = *it;
		int gains_ind = it - gains.begin();
		gains_indexes[i] = gains_ind;
		if(i > 0)
			gains_got[i] = gains_got[i] + gains_got[i-1];
		for(int j = 0; j < n; ++j)
			if(communitiesNew[gains_ind] == communitiesNew[j])
				gains[j] += 4 * Q[gains_ind][j];
			else
				gains[j] -= 4 * Q[gains_ind][j];
		communitiesNew[gains_ind] = !communitiesNew[gains_ind];
		gains[gains_ind] = gains[gains_ind] - 2*n;
	}
	vector<double>::iterator it = max_element(gains_got.begin(), gains_got.end());
	double mod_gain = *it;
	int stepsToGetMaxGain = it - gains_got.begin() + 1;
	if(mod_gain > 0)
	{
		communitiesNew = communitiesOld;
		for(int i = 0; i < stepsToGetMaxGain; ++i)
			communitiesNew[gains_indexes[i]] = !communitiesNew[gains_indexes[i]];
	}
	else
	{
		communitiesNew = communitiesOld;
		mod_gain = 0;
	}
	return mod_gain;
}

double ComboAlgorithm::Split(vector< vector<double> >& Q, const vector<double>& correctionVector, vector<int>& splitCommunity) //try to split the subnetwork with respect to the correction vector
{
	double mod_gain = 0.0;
	vector<double> sumQ = Sum(Q);
	int n = Q.size();
	for(int i = 0; i < n; ++i)
		Q[i][i] += 2 * correctionVector[i] - sumQ[i]; //adjust the submatrix
	int tries;
	if(m_num_split_attempts > 0)
		tries = m_num_split_attempts;
	else
		tries = pow(abs(log(m_current_best_gain)), m_autoC2) / m_autoC1 + 3;
	for (int tryI = 1; tryI <= tries; ++tryI)
	{
		//perform an initial simple split
		vector<int> communities(n);
		if(m_fixed_split_step > 0 && tryI <= 6 * m_fixed_split_step && tryI % m_fixed_split_step == 0) {
			//perform one of predifined split types
			int fixed_split_type = tryI / m_fixed_split_step;
			if (fixed_split_type == 1 || fixed_split_type == 2)
				communities.assign(n, 2 - fixed_split_type);
			else {
				vector<double> sum_pos = SumPos(Q, Positive);
				int node_ind;
				if (fixed_split_type == 3 || fixed_split_type == 4)
					node_ind = max_element(sum_pos.begin(), sum_pos.end()) - sum_pos.begin();
				else
					node_ind = min_element(sum_pos.begin(), sum_pos.end()) - sum_pos.begin();  
				communities.assign(n, -1);
				int community = 1;
				communities[node_ind] = community;
				while (true) {
					int next_node_ind = -1;
					double cur_min = 1e300;
					double cur_max = -1e300;
					for (int i = 0; i < n; ++i) {
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
					if (next_node_ind == -1)
						break;
					node_ind = next_node_ind;
					community ^= 1;
					communities[node_ind] = community;
				}
			}
		} else {
			for(int i = 0; i < n; ++i)
				communities[i] = int(m_bernoulli_distribution(m_random_number_generator));
		}

		double mod_gain_total = ModGain(Q, correctionVector, communities);
		double mod_gain_from_shift = 1;
		while(mod_gain_from_shift > THRESHOLD)
		{
			vector<int> communities_shifted(n);
			mod_gain_from_shift = PerformKernighansShift(Q, correctionVector, communities, communities_shifted);
			if(mod_gain_from_shift > THRESHOLD)
			{
				mod_gain_total += mod_gain_from_shift;
				communities = communities_shifted;
				if(m_debug_verify)
				{
					double mod_gain_verify = ModGain(Q, correctionVector, communities);
					if(fabs(mod_gain_verify - mod_gain_total) > THRESHOLD)
						cerr << "ERROR" << endl;
				}

			}
		}
		if(mod_gain < mod_gain_total)
		{
			splitCommunity = communities;
			mod_gain = mod_gain_total;
		}
		if(mod_gain <= 1e-6)
			tries = int(tries / 2);
	}

	if(fabs(mod_gain) < THRESHOLD)
		splitCommunity.assign(n, 1);

	return mod_gain;
}

void ComboAlgorithm::reCalc(Graph& G, vector< vector<double> >& moves, vector< vector<int> >& splits_communities, int origin, int dest)
{
	moves[origin][dest] = 0;
	if(origin != dest)
	{
		vector<int> origCommInd = G.CommunityIndices(origin);
		if(!origCommInd.empty())
		{
			vector<double> correctionVector = G.GetCorrectionVector(origCommInd, G.CommunityIndices(dest));
			vector<int> splitComunity(origCommInd.size());
			vector< vector<double> > Q = G.GetModularitySubmatrix(origCommInd);
			moves[origin][dest] = Split(Q, correctionVector, splitComunity);
			for(int i = 0; i < splitComunity.size(); ++i)
				splits_communities[dest][origCommInd[i]] = splitComunity[i];
		}
	}
}

double BestGain(const vector< vector<double> >& moves, int& origin, int& dest)
{
	double bestGain = -1;
	for(int i = 0; i < moves.size(); ++i)
		for(int j = 0; j < moves.size(); ++ j)
			if(bestGain < moves[i][j])
			{
				bestGain = moves[i][j];
				origin = i;
				dest = j;
			}
	return bestGain;
}

void DeleteEmptyCommunities(Graph& G, vector< vector<double> >& moves, vector< vector<int> >& splits_communities, int origin)
{
	if(G.DeleteCommunityIfEmpty(origin))
	{
		int commNumber = G.CommunityNumber();
		for(int i = origin; i < commNumber; ++i)
			moves[i] = moves[i+1];
		moves[commNumber].assign(commNumber+2, 0);
		for(int i = 0; i < moves.size(); ++i)
		{
			for(int j = origin; j < commNumber+1; ++j)
				moves[i][j] = moves[i][j+1];
			moves[i][commNumber+1] = 0;
		}
		for(int i = origin; i < commNumber+1; ++i)
			splits_communities[i] = splits_communities[i+1];
	}
}

void ComboAlgorithm::Run(Graph& G, int max_comunities)
{
	G.CalcModMatrix();
	G.SetCommunities(vector<int>(G.Size(), 0));
	double currentMod = G.Modularity();
	vector< vector<double> > moves(2, vector<double>(2, 0)); //results of splitting communities
	//vectors of boolean meaning that corresponding vertex should be moved to dest
	vector< vector<int> > splits_communities(2, vector<int>(G.Size(), 0)); //best split vectors
	m_current_best_gain = 1;
	int origin, dest;
	for(origin = 0; origin < G.CommunityNumber(); ++ origin)
		for(dest = 0; dest < G.CommunityNumber() + (G.CommunityNumber() < max_comunities); ++dest)
			reCalc(G, moves, splits_communities, origin, dest);

	m_current_best_gain = BestGain(moves, origin, dest);

	while(m_current_best_gain > THRESHOLD)
	{
		bool comunityAdded = dest >= G.CommunityNumber();
		G.PerformSplit(origin, dest, splits_communities[dest]);
		if(m_debug_verify)
		{
			double oldMod = currentMod;
			currentMod = G.Modularity();
			if(fabs(currentMod - oldMod - m_current_best_gain) > THRESHOLD)
				cerr << "ERROR" << endl;
		}
		if(comunityAdded && dest < max_comunities - 1)
		{
			if(dest >= moves.size() - 1)
			{
				for(int i = 0; i < moves.size(); ++i)
					moves[i].push_back(0);
				moves.push_back(vector<double>(moves.size() + 1, 0));
				splits_communities.push_back(vector<int>(G.Size(), 0));
			}
			for(int i = 0; i < dest; ++i)
			{
				moves[i][dest+1] = moves[i][dest];
				splits_communities[dest+1] = splits_communities[dest];
			}
		}

		for(int i = 0; i < G.CommunityNumber() + (G.CommunityNumber() < max_comunities); ++i)
		{
			reCalc(G, moves, splits_communities, origin, i);
			reCalc(G, moves, splits_communities, dest, i);
			if(i != dest && i < G.CommunityNumber())
				reCalc(G, moves, splits_communities, i, origin);
			if(i != origin && i < G.CommunityNumber())
				reCalc(G, moves, splits_communities, i, dest);
		}
		DeleteEmptyCommunities(G, moves, splits_communities, origin); //remove origin community if empty
		m_current_best_gain = BestGain(moves, origin, dest);
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

ComboAlgorithm::ComboAlgorithm(long long random_seed, int num_split_attempts, int fixed_split_step) :
	m_random_number_generator(random_seed),
	m_bernoulli_distribution(0.5),
	m_fixed_split_step(fixed_split_step)
{
	SetNumberOfSplitAttempts(num_split_attempts);
}

ComboAlgorithm::ComboAlgorithm(): 
	ComboAlgorithm(std::chrono::duration_cast<std::chrono::milliseconds>(
		std::chrono::steady_clock::now().time_since_epoch()).count())
{}

ComboAlgorithm::ComboAlgorithm(int num_split_attempts, int fixed_split_step) :
	ComboAlgorithm(std::chrono::duration_cast<std::chrono::milliseconds>(
		std::chrono::steady_clock::now().time_since_epoch()).count(), num_split_attempts, fixed_split_step)
{}
