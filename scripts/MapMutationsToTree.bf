
codonTo3 = {};

codonOffset = 0;

codonTo3[0] = "F";
codonTo3[1] = "L";
codonTo3[2] = "I";
codonTo3[3] = "M";
codonTo3[4] = "V";
codonTo3[5] = "S";
codonTo3[6] = "P";
codonTo3[7] = "T";
codonTo3[8] = "A";
codonTo3[9] = "Y";
codonTo3[10] = "Stop";
codonTo3[11] = "H";
codonTo3[12] = "Q";
codonTo3[13] = "N";
codonTo3[14] = "K";
codonTo3[15] = "D";
codonTo3[16] = "E";
codonTo3[17] = "C";
codonTo3[18] = "W";
codonTo3[19] = "R";
codonTo3[20] = "G";

nucCharacters = "ACGT";


/* ____________________________________________________________________	*/
/*  Converts integer identifier of codon to three-letter string of		*/
/*	nucleotides.														*/

function codeToLetters (codonCode)
{
	return nucCharacters[codonCode$16]+nucCharacters[(codonCode%16)$4]+nucCharacters[codonCode%4];	
}


/* ____________________________________________________________________	*/
function importLikelihoodFunction (dummy)
{
	SetDialogPrompt ("Choose a file containing an exported likelihood function:");
	fscanf(PROMPT_FOR_FILE, "Raw", importLF);
	
	ExecuteCommands(importLF);
	
	/* START 20070926SLKP: automatically generate the list of current LFs */
	/*	old code follows */
	/*
		fprintf (stdout, "\n\nWhat is the name (identifier) of the likelihood function in this file? : ");
		fscanf (stdin, "String", lfid);
	*/
	howManyLFs = Rows("LikelihoodFunction");
	if (howManyLFs>1)
	{
		ChoiceList  (likelihoodFnChoice,"Choose a Likelihood Function",1,NO_SKIP,LikelihoodFunction);
		if (likelihoodFnChoice<0) /* cancelled */
		{
			return 0;
		} 
	}		
	else
	{
		if (howManyLFs == 1) /* only one available - use it */
		{
			likelihoodFnChoice = 0;
		}
		else /* no LFs available; die */
		{
			fprintf (stdout, "ERROR: no likelihood functions available for ancestral state reconstruction!\n");
			return 0;
		}
	}	

	GetString (lfid,LikelihoodFunction,likelihoodFnChoice);
	/* END 20070926SLKP */

	
								/* extract identifiers from likelihood function */
	ExecuteCommands("GetString (recep, "+lfid+", -1)");
	
	dataid = (recep["Datafilters"])[0][0];		/* identifier for data filter */
	fprintf(stdout, "Retrieved data partition ", dataid, "\n");
	
	ExecuteCommands ("DataSetFilter filteredData = CreateFilter ("+dataid+",3,\"\",\"\",GeneticCodeExclusions);");

	treeid = (recep["Trees"])[0][0];				/* identifier for tree */
	fprintf(stdout, "Retrieved tree ", treeid, "\n");
	
	/* START 20070926SLKP: return 1 for success (other wise 0 */
	return 1;
	/* END 20070926SLKP */
}



/* ____________________________________________________________________	*/
/*  Assembles a key as two column vectors mapping sequence names to 	*/
/*	branch IDs.															*/

function setupMapToTree (sampling_option)
{
	/* load universal genetic code */
	skipCodeSelectionStep = 1;
	dummy = HYPHY_BASE_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"TemplateModels"+DIRECTORY_SEPARATOR+"chooseGeneticCode.def";
	ExecuteCommands ("#include \""+dummy+"\";");
	ApplyGeneticCodeTable (0);
	
	/* reconstruct ancestors from imported LF */
	if (sampling_option == 1)
	{
		ExecuteCommands("DataSet ancestralSeqs = SampleAncestors ("+lfid+");");
	} else {
		ExecuteCommands("DataSet ancestralSeqs  = ReconstructAncestors ("+lfid+");");
	}
	
	/* generate data set filter of amino acid sequences from reconstructed ancestors */
	DataSetFilter	filteredDataA  = CreateFilter(ancestralSeqs,3,"","",GeneticCodeExclusions);
	
	/* get codon frequencies from data */
	ExecuteCommands ("HarvestFrequencies (observedCEFV, "+dataid+",3,3,0);");
	HarvestFrequencies (observedCEFV, filteredData, 3, 3, 0);
	
	stateCharCount = 64;

	for (h=0; h<64; h=h+1)		/* skip stop codons */
	{
		if (_Genetic_Code[h] == 10)
		{
			stateCharCount = stateCharCount -1;
		}
	}
	ambChoice = 1;
	
	seqToBranchMap = {stateCharCount,1};	/* initialize vector */

	hShift = 0;

	_GC_ = {stateCharCount,1};
	correctCode = {stateCharCount,1};

	for (k=0; k<64; k=k+1)
	{
		if (_Genetic_Code[k]==10)
		{
			hShift = hShift+1;
		}
		else
		{
			seqToBranchMap[k-hShift] = observedCEFV[k];	/* transfer observed codon frequencies */
			_GC_[k-hShift] = _Genetic_Code[k];
			correctCode[k-hShift] = k;
		}
	}
	
	observedCEFV = seqToBranchMap;	/* reset codon frequency vector to skip termination codons */
	
	ExecuteCommands ("branchNames = BranchName ("+treeid+",-1);");
	
				/* -1 means retrieve for all branchces in tree */
	
	h = Columns (branchNames);	/* total number of branches in tree */
	
	seqToBranchMap 	= {h, 2};		/* re-using container! */
	
	/* maps sequence names to branch order in column 1 
	   and the other way around in column 2 */

	for (k=0; k<filteredData.species; k=k+1)	/* for every extant sequence (tip) */
	{
		GetString (seqName, filteredData, k);
		seqToBranchMap[k][0] = -1;
		for (v=0; v<h; v=v+1)
		{
			if (branchNames[v] % seqName)	/* string comparison */
			{
				seqToBranchMap[k][0] = v;	/* map tip name ordering to data filter */
				seqToBranchMap[v][1] = k;	/* and vice versa */
				break;
			}
		}
	}
	
	seqToBranchMap[filteredData.species][0] = h-1;
	seqToBranchMap[h-1][1] = filteredData.species;
	
	
	for (k=1; k<filteredDataA.species; k=k+1)	/* for every internal branch */
	{
		GetString (seqName, filteredDataA, k);
		seqToBranchMap[filteredData.species+k][0] = -1;
		for (v=0; v<h; v=v+1)
		{
			if (branchNames[v] % seqName)
			{
				seqToBranchMap[k+filteredData.species][0] = v;
				seqToBranchMap[v][1] = k+filteredData.species;
				break;
			}
		}
	}

	GetDataInfo    (dupInfo, filteredData);
	GetDataInfo	   (dupInfoA, filteredDataA);

	matrixTrick  = {1,stateCharCount};
	matrixTrick2 = {1,stateCharCount};

	for (h=Columns(matrixTrick)-1; h>=0; h=h-1)
	{
		matrixTrick  [h] = h;
		matrixTrick2 [h] = 1;
	}

	codonInfo  = {filteredData.species, filteredData.unique_sites};		/* nucleotide data */
	codonInfo2 = {filteredDataA.species, filteredDataA.unique_sites};	/* amino acid data */

	GetDataInfo    (dupInfo, filteredData);
	GetDataInfo	   (dupInfoA, filteredDataA);

	matrixTrick  = {1,stateCharCount};
	matrixTrick2 = {1,stateCharCount};

	for (h=Columns(matrixTrick)-1; h>=0; h=h-1)
	{
		matrixTrick  [h] = h;
		matrixTrick2 [h] = 1;
	}

	for (v=0; v<filteredData.unique_sites;v=v+1)
	{
		for (h=0; h<filteredData.species;h=h+1)
		{
			GetDataInfo (siteInfo, filteredData, h, v);
			_SITE_ES_COUNT = matrixTrick2 * siteInfo; 
			if (_SITE_ES_COUNT[0] == 1)
			{
				siteInfo = matrixTrick * siteInfo;
				codonInfo[h][v] = siteInfo[0];
			}
			else
			{
				codonInfo[h][v] = -1;
			}
		}
	}

	for (v=0; v<filteredDataA.unique_sites;v=v+1)
	{
		for (h=0; h<filteredDataA.species;h=h+1)
		{
			GetDataInfo (siteInfo, filteredDataA, h, v);
			siteInfo = matrixTrick * siteInfo;
			codonInfo2[h][v] = siteInfo[0];
		}
	}
	
	ExecuteCommands ("flatTreeRep = Abs ("+treeid+");");
	GetInformation (seqStrings, filteredData);
	
	return 0;
}




/* ____________________________________________________________________	*/
/*	Generate a tab-separated table where each row corresponds to a 		*/
/*	branch in the tree, and	each column corresponds to a codon position */
/*	in the aligned sequences.  Populate with 0's and 1's, i.e. binary 	*/
/*	variables, for presence or absence of a nonsynonymous substitution	*/
/*	at that position in each branch.									*/

function reconstructAncestors (rep, makeLabels)
{
	/* branchLengths = BranchLength (bltree, -1); */
	
	/*
	if (output_option == 0 && rep > 0)
	{
		fprintf (outfile, "\n");
	}
	*/
	
	/* START 20070926SLKP:
			Using the KEEP_OPEN option of a file does not write 
			every operation to the file immediately and close it 
			after a fprintf; this GREATLY speeds up file operations
			if many writes are performed 
			
			writeToFile is introduced as the write-to-path for 
			convenience */
	if (rep >= 0)
	{
		writeToFile = outfile+rep;
	}
	else
	{
		writeToFile = outfile;
	}
	fprintf (stdout, "Writing to file ", writeToFile, "\n");
	fprintf (writeToFile, CLEAR_FILE, KEEP_OPEN);
	
	/* END 20070926SLKP */

	if (makeLabels & output_option==0) {
	    labels = {};
	    k = filteredData.species+1;
	    p1 = seqToBranchMap[k][0];
	    pid = flatTreeRep[p1];
	    p2 = seqToBranchMap[pid][1] - filteredData.species;

	    for (site = 0; site < filteredDataA.sites; site = site+1) {
	        c1 = dupInfoA[site];
	        cd2 = codonInfo2[p2][c1];  // ancestral codon state on branch
	        aa = _GC_[cd2];
	        aa = codonTo3[aa];  // amino acid

	        labels[Abs(labels)] = aa;
	    }

        // write header row
	    for (i = 0; i < Abs(labels); i=i+1) {
	        if (i>0) {
	            fprintf(LAST_FILE_PATH, ",");
	        }
	        fprintf(LAST_FILE_PATH, labels[i], i+1);
	    }
	    fprintf(LAST_FILE_PATH, "\n");
	}



	k = filteredData.species + 1;		/* root ancestral sequence */
	
	for (h = 1; h < filteredDataA.species; h = h + 1)	/* for every internal branch */
	{
		aa_count = 0;
		
		p1 = seqToBranchMap[k][0];	/* get position of root in tree ordering */
		pid = flatTreeRep[p1];	
		p2 = seqToBranchMap[pid][1] - filteredData.species;
		
		bn = branchNames [p1];
		fprintf (writeToFile, bn);
		
		for (site = 0; site < filteredDataA.sites; site = site + 1)		/* for every site in sequence */
		{
			c1 = dupInfoA[site];
			/* START 20070926SLKP: need to handle the case 
					 when a site has all gaps - then ancestral 
					 states will also be all gaps;
					 I added this behavious to ReconstructAncestors 
					 in May 2007 */
					 
			if (codonInfo2[1][c1] >= stateCharCount)
			{
				fprintf (writeToFile, ",0");
				continue;
			}
			/* END 20070926SLKP */

			cd1 = codonInfo2[h] [c1];	/* derived codon state on branch */
			cd2 = codonInfo2[p2][c1];	/* ancestral	"		"		 */
			
			if (cd1 >= stateCharCount || cd2 >= stateCharCount)	/* --- */
			{
				fprintf (writeToFile, ",0");	/* ignore gaps */
				continue;
			}
			
			aa1 = _GC_[cd1];
			aa1 = codonTo3[aa1];		/* get amino acid character */
			aa2 = _GC_[cd2];
			aa2 = codonTo3[aa2];
			
			
			if (aa1 != aa2)		/* nonsynonymous substitution */
			{
				if (output_option == 0)		/* CSV matrix format */
				{
					fprintf (writeToFile, ",1");
				}
				else						/* list format */
				{
					fprintf (writeToFile, rep, "\t", bn, "\t", site, "\t", aa2, "-->", aa1, "\n");
					aa_count = aa_count + 1;
				}
			} else {
				if (output_option == 0)
				{
					fprintf (writeToFile, ",0");
				}
			}
		}
		/* end loop over sites */
		
		
		if (output_option == 0)
		{
			fprintf (writeToFile, "\n");	/* end of line */
		}
		else
		{
			if (aa_count == 0)
			{
				fprintf (writeToFile, rep, "\t", bn, "\tNone\n");
			}
		}
		
		k = k + 1;		/* update parent sequence */
	}
	
	/* now do the leaves */
	for (h = 0; h < filteredData.species; h = h + 1)
	{
		aa_count = 0;
		
		p1 = seqToBranchMap[h][0];
		pid = flatTreeRep[p1];
		p2 = seqToBranchMap[pid][1]-filteredData.species;
		
		bn = branchNames [p1];
		fprintf (writeToFile, bn);
		/* fprintf ("branchIDs.out", bn, "\n"); */
		/* bl = branchLengths [p1]; */
		
		
		for (site = 0; site < filteredDataA.sites; site = site + 1)
		{
			c1 = dupInfoA[site];
			c2 = dupInfo [site];
			
			cd2 = codonInfo2[p2][c1];
			cd1 = codonInfo [h] [c2];
			
			if (cd1 >= stateCharCount || cd2 >= stateCharCount)	/* --- */
			{
				fprintf (writeToFile, ",0");
				continue;
			}
			
			if (cd1 >= 0)	/* non-ambiguous codon */
			{
				aa1 = _GC_[cd1];
				aa1 = codonTo3[aa1];		/* get amino acid character */
				aa2 = _GC_[cd2];
				aa2 = codonTo3[aa2];
				
				if (aa1 != aa2)		/* nonsynonymous substitution */
				{
					if (output_option == 0)
					{
						fprintf (writeToFile, ",1");
					}
					else
					{
						fprintf (writeToFile, rep, "\t", bn, "\t", site, "\t", aa2, "-->", aa1, "\n");
						aa_count = aa_count + 1;
					}
				} 
				else 
				{
					if (output_option == 0)
					{
						fprintf (writeToFile, ",0");
					}
				}
			}
			else	/* codon contains as mixture, ignore! */
			{
				if (output_option == 0)
				{
					fprintf (writeToFile, ",0");
				}
			}
		}
		
		if (output_option == 0)
		{
			fprintf (writeToFile, "\n");	/* end line */
		}
		else
		{
			if (aa_count == 0)
			{
				fprintf (writeToFile, rep, "\t", bn, "\tNone\n");
			}
		}
	}
	/* START 20070926SLKP: CLOSE_FILE is needed to flush all to disk */
			
	fprintf (writeToFile, CLOSE_FILE);
	
	/* END 20070926SLKP */

	return 0;		/* unused */
}




/* ====================================================================================	*/
/*		MAIN LOOP																		*/
/* ==================================================================================== */


/* using universal genetic code by default */
skipCodeSelectionStep = 1;		
dummy = HYPHY_BASE_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"TemplateModels"+DIRECTORY_SEPARATOR+"chooseGeneticCode.def";
ExecuteCommands ("#include \""+dummy+"\";");
ApplyGeneticCodeTable (0);

stateCharCount = 64;	/* there are 64 possible codons, including termination codons */

for (h=0; h<64; h=h+1)	/* for every codon number */
{
	if (_Genetic_Code[h] == 10)	/* genetic code specifies stop codon (e.g. TAG) */
	{
		stateCharCount = stateCharCount -1;	/* adjust total number of non-stop codons */
	}
}


/*START 20070926SLKP: handle the case of failing to select an LF */

flag = importLikelihoodFunction(0);  // returns 1 on success, 0 otherwise
if (flag == 0) 
{
	return 0;
}
/*END 20070926SLKP*/	




ChoiceList (option, "Ancestor reconstruction option: ", 1, SKIP_NONE, 
				"Maximum likelihood", "Map most likely reconstruction to internal nodes.",
				"Sampling", "Map reconstructions at nodes in propoprtion to likelihood.");

ChoiceList (output_option, "Output option: ", 1, SKIP_NONE,
				"Binary matrix", "Write comma-separated (CSV) binary matrices to file, where 1 indicates a site-specific substitution event.",
				"List", "List amino acid substitutions per branch.");

header_option = 1;
if (output_option == 0) {
    ChoiceList(header_option, "Header option: ", 1, SKIP_NONE,
                "No", "Skip header row.",
                "Yes", "Generate column labels from amino acid sequence reconstruction at root.");
}
	
SetDialogPrompt ("Provde a filename to write output(s) to: ");
fprintf (PROMPT_FOR_FILE, CLEAR_FILE);
outfile = LAST_FILE_PATH;

/*
SetDialogPrompt ("Select a tree to extract branch lengths from: ");
fscanf(PROMPT_FOR_FILE, "Tree", bltree);
 */
 
if (option == 1)
{
	fprintf (stdout, "Enter number of replicates to sample: ");
	fscanf (stdin, "Number", max_reps);
	
	for (this_rep = 0; this_rep < max_reps; this_rep = this_rep+1)
	{
		setupMapToTree(1);
		reconstructAncestors(this_rep, header_option);
	}
}
else
{
	setupMapToTree(0);
	reconstructAncestors(-1, header_option);
}









