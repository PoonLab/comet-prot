# Bayesian graphical analysis #
## Obtaining *HyPhy* and scripts ##
*HyPhy* binaries can be downloaded for free [here](http://hyphy.org).
If you want to compile the software package, the source code can be 
obtained [here](http://github.com/veg/hyphy).
Alternatively, a POSIX-threaded *HyPhy* binary (hyphymp) for Linux 
can be obtained with the Ubuntu package manager with the command 
```
sudo apt install hyphy-pt
```
The scripts and data used here are available [here](http://github.com/PoonLab/comet-prot).
If you are running *HyPhy* from the command line, then all commands should include 
specify the path to your local installation, *e.g.*:
```
HYPHYMP BASEPATH=/usr/local/lib/hyphy <path to script>
```
Otherwise, you can run the scripts through the graphical user interface by opening 
the file through the file selection dialog. 
Example files can be found [here](https://github.com/PoonLab/comet-prot/tree/master/data)
and include the following:

```
HCV1b-NS5b.fasta
```
GenBank accession numbers retrieved from the [euHCVdb database](https://euhcvdb.ibcp.fr/euHCVdb/) for hepatitis C virus (HCV) subtype 1b nucleotide sequences with at least partial coverage of the gene NS5b.

```
HCV1b-NS5b.aliview.fa
```
Multiple sequence alignment file using MAFFT v7.305b and manually inspected and adjusted the resulting alignment with AliView v1.19-beta-3.

```
HCV1b-NS5b.mafft.fa
```
The remaining sequences in this data set ($n=536$), which ranged from 1043 to 1776 in nucleotide length.

```
HCV1b-NS5b.cleaned.fa
```
We used the built-in method in *HyPhy* to clean sequence names and remove stop codons

```
HCV1b-NS5b.phyml.nwk
```
File including the maximum likelihood phylogenetic tree was reconstructed using PhyML with a bootstrap analysis.

```
HCV1b-NS5b.lf
```
File containing the PhyML tree reconstruction.

```
HCV1b-NS5b.csv
```
Ancestral reconstruction CSV file.



## Fit codon model ##

The first step in our analysis pipeline is to fit a codon substitution model 
to the sequence alignment by running the script (which works just on command line):
```
fit_codon_model.bf
```
The following steps will be prompted:

1. Choose a genetic code

There is a large selection of genetic codes available, it will determine how nucleotide substitutions 
are interpreted as missense, nonsense or silent mutations.

2. Specify a codon data file

Enter an absolute path to the file containing the cleaned sequence alignment. It is critical that this 
alignment comprises a single open reading frame. Sequence names cannot contain any characters other 
than the alphanumeric characters and the underscore character.

3. Model options

This option determines how the model parameters are distributed across branches in the tree.
The *Local* option assigns an instance of each parameter to every branch in the tree.
The *Global* option means that each model parameter is estimated using the information from all branches in the tree.
The *Global w/variation* option models this rate variation using one of many parametric distributions, 
such as the prototypical gamma distribution.
The *Global w/variation+HM* option uses a Hidden Markov model to smooth the assignment of rate categories 
along the length of the sequence alignment.

4. Nucleotide model

The codon substitution model implemented in these scripts has a nested model of nucleotide substitution
that needs to be specified by the user (uses the 6-digit PAUP*-style model specification string).

5. Specify a tree file

Enter a relative or absolute path to the file containing the reconstructed phylogeny in a Newick tree string format
and branch length information. Do not include bootsrap information.

6. Fit a likelihood function

Specify a relative or absolute path to a file to write a serialized likelihood function, which encodes the data, model and parameter estimates. This output is written in a NEXUS format.

## Map substitutions to the tree ##

The next step in our pipeline is to reconstruct ancestral sequences in the tree based on the maximum likelihood parameter estimates of the model. On the command line, run the script:

```
MapMutationsToTree.bf
```

The following steps will be prompted:

1. Input data file

Specify a relative or absolute path file containing the likelihood function that was produced by the former script.

2. Select reconstruction option

Decide whether to sample ancestral sequences from the posterior distributions at each node of the tree.
*Sampling* enables us to accommodate the uncertainty in reconstructing ancestral states, which is exacerbated 
for ancestral nodes that are further back in time relative to the observed sequences.

3. Output options

The first option is to generate a binary matrix where each row corresponds to a branch of the tree, and each 
column corresponds to a codon site in the alignment. This matrix is written to the output file in a comma-separated 
tabular (CSV) format. Then the user is prompted to specify if they want the CSV to begin with a header row.
The second option is to output a tab-delimited tabular file where each row corresponds to a inferred non-synonymous substitution, and the columns correspond to the replicate number (if sampling), branch label, site, and the ancestral and derived amino acids.

## BGM analysis ##

Perform a Bayesian graphical model analysis and works just on command line:

```
bayesgraph.ib
```
which depends on the utility:
```
bayesgraph.ibf
```

The following steps will be prompted:

1. Input data matrix

Specify a relative or absolute path to the CSV file containing the substitution map matrix that was produced by the former script.

2. Header options

Specify wether your input file has a header.

3. Filter sites

Specify the minimum number of substitutions for a site to be included in the BGM model.

4. MCMC settings

There are four settings that the user needs to specify for running a Markov chain Monte Carlo (MCMC) sample.
First, the user has to specify the maximum number of parents that will be allowed per node.
Second, the user needs to indicate the number of steps to discard as a *burn-in* period.
Third, the user needs to specify the number of steps to run the chain sample following the *burn-in* period.
Lastly, the user must specify the number of steps to extract from this chain sample.

The script produces three kinds of outputs.
First, the script will output the marginal posterior probabilities for directed edges as a CSV formatted file.
Next, the script will write this consensus BGM using the network visualization language DOT, 
which can be converted into an image by several programs such as GraphViz, Cytoscape and Gephi.
Finally, the script will record the posterior probability trace for all steps sub-sampled from the original MCMC sample.











