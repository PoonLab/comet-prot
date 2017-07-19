/*_____________________________________________________________________ */
function computeScalingFactorB(rateMatrix, baseFreqs)
{
	B = 0;
	for (n1 = 0; n1 < Rows(rateMatrix); n1 = n1+1)
	{
		for (n2 = 0; n2 < Columns(rateMatrix); n2 = n2+1)
		{
			if (n2 != n1)
			{
				B = B + baseFreqs[n1]*baseFreqs[n2]*rateMatrix[n1][n2];
			}
		}
	}
	return B;
}


/*_____________________________________________________________________ */

NICETY_LEVEL = 3;
VERBOSITY_LEVEL = 0;

dummy = HYPHY_BASE_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR + "TemplateModels" + DIRECTORY_SEPARATOR + "chooseGeneticCode.def";
ExecuteCommands ("#include \""+dummy+"\";");


SetDialogPrompt ("Please specify a codon data file:");
DataSet codon_ds = ReadDataFile (PROMPT_FOR_FILE);
DataSetFilter codon_dsf = CreateFilter (codon_ds,3,"","",GeneticCodeExclusions);
fprintf (stdout,"\n______________READ THE FOLLOWING DATA______________\n", codon_ds);

	/* generate codon model (MG94customModel) */
#include "fit_codon_model.ibf";


	/* read in tree from file */
ACCEPT_BRANCH_LENGTHS 	= 1;
ACCEPT_ROOTED_TREES		= 1;
SetDialogPrompt ("Please select a file containing a tree with branch lengths: ");
fscanf(PROMPT_FOR_FILE, "String", tree_string);

Tree	codon_tree = tree_string;



	/* constrain branch lengths */
global scalingB 		= computeScalingFactorB (MG94custom, vectorOfFrequencies);

branchNames 	= BranchName (codon_tree, -1);
branchLengths	= BranchLength (codon_tree, -1);

for (k = 0; k < Columns(branchNames)-1; k=k+1)
{
	ExecuteCommands("codon_tree." + branchNames[k] + ".synRate:=" + branchLengths[k] + "/scalingB;");
}


SetDialogPrompt ("Please specify a file to export likelihood function: ");
fprintf (PROMPT_FOR_FILE, CLEAR_FILE);
LikelihoodFunction codon_lf = (codon_dsf, codon_tree);

AUTO_PARALLELIZE_OPTIMIZE = 1;	/* attempt to use MPI */
Optimize (res, codon_lf);
AUTO_PARALLELIZE_OPTIMIZE = 0;

LIKELIHOOD_FUNCTION_OUTPUT = 7;		/* save LF to file */
fprintf(LAST_FILE_PATH, codon_lf);
LIKELIHOOD_FUNCTION_OUTPUT = 2;		/* reset to default (?) */


fprintf (stdout, "\n______________RESULTS______________\n",codon_lf);
