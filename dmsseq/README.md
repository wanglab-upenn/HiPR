
# HiPR - High-throughput Probabilistic inference of RNA secondary structures 

Estimate secondary structure and base pairing posteriors for a given RNA sequence based on the distribution of read fragments along the locus.

## Input
This program requires a file containing the sequence and initial secondary structure of an RNA of interest (see e.g., `test.starting_structures`), a file containing reverse trascriptase (RT) stop-based (DMS-Seq/Structure-seq) read-out (e.g., `test.reads`), and a file containing the initial estimates of per-nucleotide modification rates (e.g., `test.rates`). A Bayesian MCMC algorithm is then used to estimate the base pairing posterior that best fits the observed sequencing reads.

## Output
The results and intermediate files are written to a directory (HiPR_output/ by default).

The output file `HiPR_posterior.txt` contains the base pairing posteriors at each nucleotide position, one entry per line.

The output file `HiPR_structure.txt` contains the consensus secondary structure.

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

Perl (v5.x). Important: Forks::Super module is required to run HiPR. To install, e.g., perl -MCPAN -e "install 'Forks::Super'" or if it fails to install because some tests are failing during installation, these tests can be skipped altogether: perl -MCPAN -e "notest('install','Forks::Super')"


Boost C++ libraries  

Standard POSIX programs (awk, grep, bash)

C++ compiler (for re-compilation if necessary)

===========================================================================
## Description of the files included into HiPR distribution

HiPR.sh -- main HiPR pipeline

HiPR_MCMC -- Linux executable (may need to be recompiled)

HiPR_MCMC.mac -- Mac executable

HiPR_MCMC_wrapper.pl -- wrapper script 

INSTALL.md -- instructions for re-compilation and installation of required modules

README.md -- this file

dms_get_majority_structure.sh -- auxiliary script for computing consensus structure

dms_mcmc_collect_structures.sh -- auxiliary script for collecting sampled structures

test.rates -- example file with per-nucleotide DMS modification rates (ACUG)

test.reads -- example file in collapsed reads format

test.start_structures -- example file with sequence and initial structure estimate(s)

## File formats

### Rates file:
Two lines of comma-separated modification rates for A, C, U, and G nucleotides for paired (first line) and unpaired states (second line), e.g.,:
```
0.0278,0.0112,0.0064,0.0091
0.2027,0.0844,0.0110,0.0099
````
Modification rates for unpaired nucleotides are assumed to be always greater than modification rates for paired positions.
Also see `test.rates` for an example.

### Reads file:
Contains observed read intervals.
Each line (read interval) has the following format:
```
read<ReadNumber>@@<ReadCount>@@<ReadLength> <ReadStart0> <ReadEnd0> 3
```
E.g.,
```
read1@@1@@33	0	32	3
read2@@6@@37	0	36	3
read3@@1@@40	0	39	3
...
```
where the second read interval (read2) corresponds to 6 reads with length=37 that span positions `0...36`. Also see `test.reads` for an example.

### Starting structures file
Contains one or more sequences and starting structures:
```
<sequence1>
<secondary-structure-in-parenthesis-format1>
<sequence2>
<secondary-structure-in-parenthesis-format2>
...
```
E.g.,
```
UAACUCUAAUUCCCAUUUUGCAAAAUUUCCAGUACCUUUGUCACAAUCCUAACACAUUAUCGGGAGCAGUGUCUUCCAUAAUGUAUAAAGAACAAGGUAGUUUUUACCUACCACAGUGUCUGUAUCGGAGACAGUGA
..............(((.((((((............)))).)).))).((...(((((((.(((((......))))))))))))....)).....(((((.......))))).((.(((((.......))))).)). (-30.00)
UAACUCUAAUUCCCAUUUUGCAAAAUUUCCAGUACCUUUGUCACAAUCCUAACACAUUAUCGGGAGCAGUGUCUUCCAUAAUGUAUAAAGAACAAGGUAGUUUUUACCUACCACAGUGUCUGUAUCGGAGACAGUGA
..............(((.((((((............)))).)).))).((...(((((((.(((((......))))))))))))....)).....((((((....).))))).((.(((((.......))))).)). (-30.00)
```
(also see `test.starting_structures` for an example)

