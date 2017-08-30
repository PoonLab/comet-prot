VERBOSITY_LEVEL=2;

ExecuteAFile ("bayesgraph.ibf");

//data = import_data("../data/baynet.csv", 1);
//data = import_data("../data/toBGM-AD.csv", 1);
//fprintf (stdout, names, "\n");

data = import_data("../data/1a_1b.anc.csv", 1);
data = filter_data_matrix(data, 10);  // drop sites with no substitutions

ncol = Columns(data);
fprintf (stdout, "Reduced data to ", ncol, " sites\n");

// assign a node for each codon site in the data matrix
nodes = {};  // initialize associative list (dictionary)

//FIXME: `names` is undefined!

// for every column
for (i=0; i < ncol; i=i+1) {
    // add_discrete_node creates an entry into assoc.list to represent a variable
    nodes[Abs(nodes)] = add_discrete_node(names[i], 1, 0, 2);  // node_id, max_parents, sample_size, nlevels
}

// construct the BGM object
BayesianGraphicalModel my_bgm = (nodes);


// assign our data matrix to the BGM object
attach_data("my_bgm", data, 0, 0, 0);


VERBOSITY_LEVEL=7;  // report everything!

// run MCMC sample (total number of steps, burnin, number of samples)
result = order_MCMC("my_bgm", 100000, 10000, 100);
//fprintf(stdout, result);

/*
write_edgelist("../results/BGM-AD.edgelist2.withSubtype.csv", result, 1);
//mcmc_graph_to_dotfile("../results/BGM.dot", 0.8, result, nodes);
write_MCMC_chain("../results/BGM-AD.chain2.withSubtype.csv", result);
*/
