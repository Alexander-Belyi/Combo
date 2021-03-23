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
#include <iostream>
#include <optional>
#include <string>
using namespace std;

int main(int argc, char** argv)
{
	optional<size_t> max_communities = nullopt;
	string file_suffix = "comm_comboC++";
	// Modularity Resolution Parameter
	// as per Newman 2016 (https://journals.aps.org/pre/abstract/10.1103/PhysRevE.94.052315)
	double mod_resolution = 1.0;
	int num_split_attempts = 0;
	int fixed_split_step = 0;
	bool start_separate = false;
	bool treat_as_modularity = false;
	int info_output_level = 0;
	optional<string> intermediate_results_file_name = nullopt;
	if (argc < 2) {
		cerr << "Error: provide path to edge list (.edgelist) or pajeck (.net) file" << endl;
		return -1;
	}
	string file_name = argv[1];
	if (argc > 2 && string(argv[2]) != "INF")
		max_communities = atoi(argv[2]);
	if (argc > 3) 
		mod_resolution = atof(argv[3]);
	if (argc > 4) 
		file_suffix = argv[4];
	if (argc > 5)
		num_split_attempts = atoi(argv[5]);
	if (argc > 6)
		fixed_split_step = atoi(argv[6]);
	if (argc > 7)
		start_separate = atoi(argv[7]);
	if (argc > 8)
		treat_as_modularity = atoi(argv[8]);
	if (argc > 9)
		info_output_level = atoi(argv[9]);
	if (argc > 10)
		intermediate_results_file_name = argv[10];
	Graph graph = ReadGraphFromFile(file_name, mod_resolution, treat_as_modularity);
	if (graph.Size() <= 0) {
		cerr << "Error: graph is empty" << endl;
		return -1;
	}
	ComboAlgorithm combo(nullopt, num_split_attempts, fixed_split_step, info_output_level);
	combo.Run(graph, max_communities, start_separate, intermediate_results_file_name);
	graph.PrintCommunity(file_name.substr(0, file_name.rfind('.')) + "_" + file_suffix + ".txt");
	cout << graph.Modularity() << endl;
	return 0;
}
