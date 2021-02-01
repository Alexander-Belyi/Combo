# Combo

## Description

This is an implementation (for Modularity maximization) of the community
detection algorithm called "Combo" described in the paper "General optimization
technique for high-quality community detection in complex networks" by
Stanislav Sobolevsky, Riccardo Campari, Alexander Belyi and Carlo Ratti.
If you have any feedback, bug reports or questions, please submit an issue or start a discussion.

If you use this code, please, consider citing:

```
Sobolevsky S., Campari R., Belyi A., and Ratti C. "General optimization technique for high-quality community detection in complex networks" Phys. Rev. E 90, 012811
```

## Running

First you'll need to compile the program. If you're on a Linux/Mac box, make
sure you have `g++` installed, and then run

```
make
```

Once compiled, you can run the program with (parameters' order matters)

```
./comboCPP path_to_network_file.net [max_number_of_communities] [mod_resolution] [file_suffix] [num_split_attempts] [fixed_split_step]
```

The options are as follows:
* `path_to_network_file` - path to the file in Pajek .net format
* `max_number_of_communities` - maximal number of communities to be found
  (default is "INF" for infinite)
* `mod_resolution` - modularity resolution parameter (default is 1)
* `file_suffix` - suffix appended to output file (default is "comm_comboC++")
* `num_split_attempts` - number of attempts to split each community on every iteration (default is 0 meaning this number will be selected automatically on each iteration)
* `fixed_split_step` - number of steps between trying to apply each of 6 predifined splits (default is 1 - first 6 split attempts will be of predefined types, use 0 to use only random splits)

For example, you can make sure the compilation worked correctly by running:
```
./comboCPP karate.net
```
on the included `karate.net` file.

## The Output

The program outputs one file named `path_to_network_file_<file_suffix>.txt`
containing numbers of communities for each vertex on separete line.  And it
writes achieved modularity score to standart output.
