VERBOSITY_LEVEL=2;

// import helper functions from another file
ExecuteAFile ("bayesgraph.ibf");

SetDialogPrompt("Select file containing substitution matrix in CSV format");
fprintf(PROMPT_FOR_FILE, KEEP_OPEN);

fprintf(stdout, "Does this file contain a header row? (Y\\n)");
fscanf(stdin, "String", has_header);
header_option = 1;
if (has_header != 0 && (has_header=="n" || has_header=="N")) {
    header_option = 0;  // will generate arbitrary labels
}
fprintf(stdout, "header_option: ", header_option, "\n");

data = import_data(LAST_FILE_PATH, header_option);
fprintf(LAST_FILE_PATH, CLOSE_FILE);


/* determine the maximum number of 1's per column */
max_count = 0;
for (i=0 ; i < Columns(data); i=i+1) {
    count = 0;
    for (j=0; j<Rows(data); j=j+1) {
        if (data[j][i] < 0 || data[j][i] > 1) {
            fprintf (stdout, "ERROR: non-binary value detected at (", j, ",", i, ")\n");
            fprintf (stdout, "If you want to analyze factors with more than two levels, then you need to customize this script at the 'add_discrete_node' step and comment out this test.\n");
            return 0;
        }
        count = count + data[j][i];
    }
    if (count > max_count) {
        max_count = count;
    }
}






while(1) {
    fprintf(stdout, "Enter minimum number of substitutions per site (1-", max_count, "): ");
    fscanf(stdin, "Number", cutoff);
    if (cutoff > 0 && cutoff <= max_count) {
        break;
    }
    fprintf(stdout, "Cutoff must be greater than 0 and less than or equal to ", max_count, "\n");
}

data = filter_data_matrix(data, cutoff);  // drop sites with no substitutions

ncol = Columns(data);
if (ncol<=1) {
    fprintf (stdout, "Oops, your cutoff reduced the number of columns to one!\n");
    return 0;
}
fprintf (stdout, "Reduced data to ", ncol, " sites\n");


fprintf(stdout, "Maximum number of parents per node (default 1):");
fscanf(stdin, "Number", max_parents);
if (max_parents <= 0) {
    max_parents = 1;
}
if (max_parents > 3) {
    fprintf (stdout, "***********************************************************\n");
    fprintf (stdout, " WARNING: You have requested more than 3 parents per node. \n");
    fprintf (stdout, " The analysis will run, but you should be aware that model \n");
    fprintf (stdout, " complexity increases super-exponentially with this param- \n");
    fprintf (stdout, " eter and you should be skeptical of the outputs, especia- \n");
    fprintf (stdout, " lly if the chain sample is too short.                     \n");
    fprintf (stdout, "***********************************************************\n");
}

/* configure BGM analysis */
fprintf(stdout, "Length of burn-in (default 10000): ");
fscanf(stdin, "Number", burnin);
if (burnin <= 0) {
    burnin = 10000;
}
fprintf(stdout, "Length of chain post-burn-in (default 100000): ");
fscanf(stdin, "Number", chain);
if (chain <= 0) {
    chain = 100000;
}
fprintf(stdout, "Number of samples (default 100): ");
fscanf (stdin, "Number", nsamples);
if (nsamples <= 0) {
    nsamples = 100;
}
if (nsamples > chain) {
    fprintf (stdout, "Warning: requested more samples than length of chain.\n");
    nsamples = chain;
}


SetDialogPrompt("Provide a file to write edge posterior outputs");
fprintf(PROMPT_FOR_FILE, CLEAR_FILE);
dot_file = LAST_FILE_PATH + ".dot";
trace_file = LAST_FILE_PATH + ".trace.csv";

// assign a node for each codon site in the data matrix
nodes = {};  // initialize associative list (dictionary)


// for every column
for (i=0; i < ncol; i=i+1) {
    // add_discrete_node creates an entry into assoc.list to represent a variable
    nodes[Abs(nodes)] = add_discrete_node(names[i], max_parents, 0, 2);  // node_id, max_parents, sample_size, nlevels
}

// construct the BGM object
BayesianGraphicalModel my_bgm = (nodes);


// assign our data matrix to the BGM object
attach_data("my_bgm", data, 0, 0, 0);


VERBOSITY_LEVEL=7;  // report everything!

// run MCMC sample (total number of steps, burnin, number of samples)
result = order_MCMC("my_bgm", chain, burnin, nsamples);
//fprintf(stdout, result);

write_edgelist(LAST_FILE_PATH, result, Abs(nodes), 1);
mcmc_graph_to_dotfile(dot_file, 0.8, result, nodes);
write_MCMC_chain(trace_file, result, nsamples);




