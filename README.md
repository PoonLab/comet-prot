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
the file through the file selection dialog. Example files can be found [here](https://github.com/PoonLab/comet-prot/tree/master/data).
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





