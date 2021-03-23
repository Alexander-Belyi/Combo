# Combo

## Updates

1. Added new parameters to track progress.
2. Adde support for `.csv` files containing adjacency matrix.
3. Added `treat_as_modularity` flag to indicate that we should not calculate modularity matrix and use received edge weights instead. This allows users to use __Combo__ to solve general _clique partitioning problem_.
Combining this option with reading a matrix from csv allows to apply __Combo__ to any custom _modularity_ matrix.

4. We have a python wrapper around this code available here: [pyCombo](https://github.com/Casyfill/pyCombo).

## Description

This is an implementation (for Modularity maximization) of the community detection algorithm called "Combo" described in the paper "General optimization technique for high-quality community detection in complex networks" by Stanislav Sobolevsky, Riccardo Campari, Alexander Belyi and Carlo Ratti.
If you have any feedback, bug reports or questions, please submit an issue or start a discussion.

If you use this code, please, consider citing:

```
Sobolevsky, S., Campari, R., Belyi, A. and Ratti, C., 2014. General optimization technique for high-quality community detection in complex networks. Physical Review E, 90(1), p.012811.
```

## Running

First you will need to compile the program. If you're on a Linux/Mac box, make
sure you have `g++` installed, and then run

```
make
```

Once compiled, you can run the program with (parameters' order matters):

```
./comboCPP path_to_network_file.net [max_number_of_communities] [mod_resolution] [file_suffix] [num_split_attempts] [fixed_split_step] [start_separate] [treat_as_modularity] [info_output_level] [intermediate_results_file_name]
```

The options are as follows:
* `path_to_network_file` - path to the file in Pajek `.net` format, `.edgelist` file with list of edges or `.csv` file with adjacency matrix  
* `max_number_of_communities` - maximal number of communities to be found
  (default is "INF" for infinite)
* `mod_resolution` - modularity resolution parameter (default is 1)
* `file_suffix` - suffix appended to output file (default is "comm_comboC++")
* `num_split_attempts` - number of attempts to split each community on every iteration (default is 0 meaning this number will be selected automatically on each iteration)
* `fixed_split_step` - number of steps between trying to apply each of 6 predefined splits (default is 1 - first 6 split attempts will be of predefined types, use 0 to use only random splits)
* `start_separate` - 0 or 1 (default is 0), indicating if Combo should start from assigning each node into its own separate community. This can help to achieve higher modularity, but makes running time much longer for large networks
* `treat_as_modularity` - 0 or 1 (default is 0), indicating that edge weights will be treated as modularity scores
* `info_output_level` - 0 (default) or higher, indicating how much progress info should be printed out
* `intermediate_results_file_name` - path to the file where intermediate community assignment will be saved, if empty (default) no intermediate results will be saved.

For example, you can make sure the compilation worked correctly by running:
```
./comboCPP karate.net
```
on the included `karate.net` file.

## The Output

The program outputs one file named `path_to_network_file_<file_suffix>.txt` containing community labels for each vertex on a separate line.  And it writes achieved modularity score to standard output.
