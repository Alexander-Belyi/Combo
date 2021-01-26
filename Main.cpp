/*
    Copyright 2014
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

#include "Graph.h"
#include "Combo.h"
#include <ctime>
#include <string>
#include <iostream>
using namespace std;

#define INF 1000000000

extern bool use_fixed_tries;

int main(int argc, char** argv)
{
	int max_comunities = INF;
	string file_suffix = "comm_comboC++";
	// Modularity Resolution Parameter
	// as per Newman 2016 (https://journals.aps.org/pre/abstract/10.1103/PhysRevE.94.052315)
	double mod_resolution = 1.0;
	if(argc < 2)
	{
		cerr << "Error: provide path to edge list (.edgelist) or pajeck (.net) file" << endl;
		return -1;
	}
	if(argc > 2)
	{
		if(string(argv[2]) != "INF")
		max_comunities = atoi(argv[2]);
	}
	if(argc > 3) 
	{
		mod_resolution = atof(argv[3]);
	}
	if(argc > 4) 
	{
		file_suffix = argv[4];
	}
	if(argc > 5)
		use_fixed_tries = atoi(argv[5]);


	string fileName = argv[1];
	srand(time(0));

	Graph G;
	string ext = fileName.substr(fileName.rfind('.'), fileName.length() - fileName.rfind('.'));
	if(ext == ".edgelist")
		G.ReadFromEdgelist(fileName, mod_resolution);
	else if(ext == ".net")
		G.ReadFromPajeck(fileName, mod_resolution);
	if(G.Size() <= 0)
	{
		cerr << "Error: graph is empty" << endl;
		return -1;
	}

	clock_t startTime = clock();
	RunCombo(G, max_comunities);

	//cout << fileName << " " << G.Modularity() << endl;
	//cout << "Elapsed time is " << (double(clock() - startTime)/CLOCKS_PER_SEC) << endl;

	G.PrintCommunity(fileName.substr(0, fileName.rfind('.')) + "_" + file_suffix + ".txt");
	cout << G.Modularity() << endl;
	return 0;
}
