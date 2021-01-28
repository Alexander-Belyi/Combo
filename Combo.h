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
#include <random>
#include <vector>

class ComboAlgorithm {
public:
    ComboAlgorithm();
    explicit ComboAlgorithm(long long random_seed, int num_split_attempts = 0, int fixed_split_step = 0);
    ComboAlgorithm(int num_split_attempts, int fixed_split_step);
    void Run(Graph& G, int max_comunities);
    void SetFixedSplitStep(int fixed_split_step) {m_fixed_split_step = fixed_split_step;}
    void SetNumberOfSplitAttempts(int split_tries);
private:
    //settings
    const bool m_debug_verify = false;
    // number split attempts; 0 - autoadjust this number based on m_current_best_gain
    int m_num_split_attempts;
    // step number to apply predifined split; 0 - use only random splits
    // if >0 sets up the usage of 6 fixed type splits on every m_fixed_split_step
    int m_fixed_split_step;
    double m_autoC1;
    double m_autoC2;
    //implementation
    std::mt19937 m_random_number_generator;
    std::bernoulli_distribution m_bernoulli_distribution;
    double m_current_best_gain;
    void reCalc(Graph& G, std::vector< std::vector<double> >& moves, std::vector< std::vector<int> >& splits_communities, int origin, int dest);
    double Split(std::vector< std::vector<double> >& Q, const std::vector<double>& correctionVector, std::vector<int>& splitCommunity);
};
