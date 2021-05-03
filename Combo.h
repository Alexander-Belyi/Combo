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

#ifndef COMBO_H
#define COMBO_H

#include "Graph.h"
#include <cstdint>
#include <optional>
#include <random>
#include <string>
#include <vector>

class ComboAlgorithm {
public:
    explicit ComboAlgorithm(std::optional<uint_fast32_t> random_seed = std::nullopt,
        int num_split_attempts = 0, int fixed_split_step = 0, int output_info_level = 0);
    void Run(Graph& graph, std::optional<size_t> max_communities = std::nullopt, bool start_separate = false,
        std::optional<std::string> intermediate_result_file_name = std::nullopt);
    void SetFixedSplitStep(int fixed_split_step) {m_fixed_split_step = fixed_split_step;}
    void SetNumberOfSplitAttempts(int split_tries);
    void SetOutputInfoLevel(int output_info_level) {m_output_info_level = output_info_level;}
private:
    //settings
    const bool m_debug_verify = false;
    // number of split attempts; 0 - autoadjust this number based on m_current_best_gain
    int m_num_split_attempts;
    // step number to apply predefined split; 0 - use only random splits
    // if >0 sets up the usage of 6 fixed type splits on every m_fixed_split_step
    int m_fixed_split_step;
    int m_output_info_level;

    //implementation
    double m_autoC1;
    double m_autoC2;
    std::mt19937 m_random_number_generator;
    std::bernoulli_distribution m_bernoulli_distribution;
    double m_current_best_gain;
    void ReCalc(Graph& graph, std::vector< std::vector<double> >& moves,
        std::vector< std::vector<bool> >& splits_communities, size_t origin, size_t destination);
    double Split(std::vector< std::vector<double> >& Q, const std::vector<double>& correction_vector, std::vector<int>& to_be_moved);
};

#endif //COMBO_H
