
# HiPR - High-throughput Probabilistic inference of RNA secondary structures 

Estimate secondary structure and base pairing posteriors for a given RNA sequence based on the distribution of read fragments along the locus.

## Input
This program requires a file containing the sequence and initial secondary structure of an RNA of interest, a file containing DMS-seq reads along the locus, and a file containing the initial estimates of per-nucleotide modification rates. A Bayesian MCMC algorithm is then used to estimate the base pairing posterior that best fits the observed sequencing reads.

## Output
The results and intermediate files are written to a directory (HiPR_output/ by default).

The output file HiPR_posterior.txt contains the base pairing posteriors at each nucleotide position, one entry per line.

The output file HiPR_structure.txt contains the consensus secondary structure.

=========================================================================
## HiPR main script (HiPR.sh)

=========================================================================
### Operating Instructions:

Documentation:
perldoc HiPR.sh
pod2text HiPR.sh

### Usage:
```
HiPR.sh reads_file rates_file structure_file [-outDir OutputDir] [-locusName LocusName] [-n Niter] [-numCPU Ncpus]  [-rmin MinReadLength] [-rmax MaxReadLength]"
```

### Example of running HiPR using included sample files:
```
bash HiPR.sh test.reads test.rates test.start_structures -rmin 15 -rmax 45 -n 1000 -locusName test
```

If the above executes successfully, concluding messages should be displayed:

```
HiPR algorithm complete

Consensus secondary structure:
HiPR_output/test.HiPR_structure.txt

Base-pairing posteriors at each nucleotide:
HiPR_output/test.HiPR_posterior.txt
```

===========================================================================
## Software Requirements:
Linux-based operating system, or Mac OS.

Perl (v5.x). Important: Forks::Super module is required to run HiPR. To install, e.g., perl -MCPAN -e "install 'Forks::Super'"

Boost C++ libraries  

Standard POSIX programs (awk, grep, bash)

C++ compiler (for re-compilation if necessary)

===========================================================================
## Description of the files included into HiPR distribution

HiPR.sh -- main HiPR pipeline

HiPR_MCMC -- Linux executable (may need to be recompiled)

HiPR_MCMC.mac -- Mac executable

HiPR_MCMC_wrapper.pl -- wrapper script 

INSTALL -- instructions for re-compilation and installion of required modules

README -- this file

dms_get_majority_structure.sh -- auxiliary script for computing consensus structure

dms_mcmc_collect_structures.sh -- auxiliary script for collecting sampled structures

test.rates -- example file with per-nucleotide DMS modification rates (ACUG)

test.reads -- example file in collapsed reads format

test.start_structures -- example file with sequence and initial structure estimate(s)
