
ExecuteAFile ("bayesgraph.ibf.txt");

//data = import_data("../data/baynet.csv", 1);
data = import_data("../data/toBGM-AD.csv", 1);
//fprintf (stdout, names, "\n");

nodes = {};
for (i=0; i < 15; i=i+1) {
    // node_id, max_parents, sample_size, nlevels
    nodes[Abs(nodes)] = add_discrete_node(names[i], 3, 0, 2);
}
nodes[Abs(nodes)] = add_discrete_node(names[Abs(nodes)], 3, 0, 3);  // region (origin)
nodes[Abs(nodes)] = add_discrete_node(names[Abs(nodes)], 4, 0, 2);  // baseline/failures
nodes[Abs(nodes)] = add_discrete_node(names[Abs(nodes)], 3, 0, 2);  // simplified SCUEAL subtype

BayesianGraphicalModel my_bgm = (nodes);

/*
constraints = {Abs(nodes), Abs(nodes)};
for (i=0; i<15; i=i+1) {
    constraints[i][15] = -1;  // banned edge
}
setConstraints("my_bgm", constraints);
*/

attach_data("my_bgm", data, 0, 0, 0);
VERBOSITY_LEVEL=7;
result = order_MCMC("my_bgm", 100000, 100000, 2000);
//fprintf(stdout, result);

write_edgelist("../results/BGM-AD.edgelist2.withSubtype.csv", result, 1);
//mcmc_graph_to_dotfile("../results/BGM.dot", 0.8, result, nodes);
write_MCMC_chain("../results/BGM-AD.chain2.withSubtype.csv", result);
