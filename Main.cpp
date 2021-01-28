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
#include <ctime>
#include <iostream>
#include <string>
using namespace std;

int main(int argc, char** argv)
{
	int max_comunities = 2e9;
	string file_suffix = "comm_comboC++";
	// Modularity Resolution Parameter
	// as per Newman 2016 (https://journals.aps.org/pre/abstract/10.1103/PhysRevE.94.052315)
	double mod_resolution = 1.0;
	int num_split_attempts = 0;
	int fixed_split_step = 0;
	if (argc < 2) {
		cerr << "Error: provide path to edge list (.edgelist) or pajeck (.net) file" << endl;
		return -1;
	}
	if (argc > 2) {
		if(string(argv[2]) != "INF")
			max_comunities = atoi(argv[2]);
	}
	if (argc > 3) 
		mod_resolution = atof(argv[3]);
	if (argc > 4) 
		file_suffix = argv[4];
	if (argc > 5)
		num_split_attempts = atoi(argv[5]);
	if (argc > 6)
		fixed_split_step = atoi(argv[5]);

	string fileName = argv[1];
	srand(time(0));

	Graph G(mod_resolution);
	string ext = fileName.substr(fileName.rfind('.'), fileName.length() - fileName.rfind('.'));
	if(ext == ".edgelist")
		G.ReadFromEdgelist(fileName);
	else if(ext == ".net")
		G.ReadFromPajeck(fileName);
	if(G.Size() <= 0)
	{
		cerr << "Error: graph is empty" << endl;
		return -1;
	}

	clock_t startTime = clock();
	ComboAlgorithm combo(num_split_attempts, fixed_split_step);
	combo.Run(G, max_comunities);

	//cout << fileName << " " << G.Modularity() << endl;
	//cout << "Elapsed time is " << (double(clock() - startTime)/CLOCKS_PER_SEC) << endl;

	G.PrintCommunity(fileName.substr(0, fileName.rfind('.')) + "_" + file_suffix + ".txt");
	cout << G.Modularity() << endl;
	return 0;
}
