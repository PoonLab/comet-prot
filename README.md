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
