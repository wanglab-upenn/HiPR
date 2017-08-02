set -e

# functions
function get_ncpus() {
   # Linux, MAC
   ncpu=$(getconf _NPROCESSORS_ONLN 2>/dev/null)
   # FreeBSD
   [ -z "$ncpu" ] && ncpu=$(getconf NPROCESSORS_ONLN)
   # Solaris
   [ -z "$ncpu" ] && ncpu=$(ksh93 -c 'getconf NPROCESSORS_ONLN')
   # Windows
   [ -z "$ncpu" ] && ncpu=$(echo %NUMBER_OF_PROCESSORS%)
   # Other... 
   [ -z "$ncpu" ] && ncpu=1  

   echo "${ncpu}"
}



function die() {
   echo "$1" >&2
   exit 1
}

VIEWDOC=pod2text
command -v "${VIEWDOC}" > /dev/null 2>&1 || { VIEWDOC=""; }

scriptDir=`dirname $0`
PATH=$PATH:${scriptDir}
mcmcCmd=HiPR_MCMC
if [ $OSTYPE == "darwin"* ]; then
  mcmcCmd=HiPR_MCMC.mac
fi
mcmcCmdWrapper=HiPR_MCMC_wrapper.pl


# check HiPR executables
HIPRMCMC=${HIPRMCMC:-$(type -P ${mcmcCmd})}
command -v "${HIPRMCMC}" > /dev/null 2>&1 || { echo >&2 "$0 requires ${mcmcCmd}, but it is not found. Please set HIPRMCMC or make sure ${mcmcCmd} is in the path"; exit 1; }
HIPRMCMCwrapper="${scriptDir}/${mcmcCmdWrapper}"
if [ ! -s "${HIPRMCMCwrapper}" ]; then 
  echo >&2 "$0 requires ${mcmcCmdWrapper}, but it is not found. Please make sure ${mcmcCmdWrapper} is in the path";
  exit 1;
fi

# required inputs [must be specified first, before optional arguments]
COLLAPSEDREADFILE=$1 # DMS.collapsed.reads
MODRATEFILE=$2 # per-nucleotide modification rates for unpaired/paired bases
STARTSTRUCTFILE=$3 # e.g., structures sampled from Boltzmann distribution

if [ $# -lt 3 ]; then
  #echo >&2 "USAGE: $0 reads_file rates_file structure_file [OutputDir] [LocusName] [Niter] [Ncpus] [ModelType] [MinReadLength] [MaxReadLength]"
  echo >&2 "USAGE: $0 reads_file rates_file structure_file [-outDir OutputDir] [-locusName LocusName] [-n Niter] [-numCPU Ncpus]  [-rmin MinReadLength] [-rmax MaxReadLength]"
  if [ -n "${VIEWDOC}" ]; then
    ${VIEWDOC} $0;
  fi
  exit 1;
fi

if [ ! -s "${COLLAPSEDREADFILE}" ]; then
  die "${COLLAPSEDREADFILE} does not exist or is empty!"
fi
if [ ! -s "${MODRATEFILE}" ]; then
  die "${MODRATEFILE} does not exist or is empty!"
fi
if [ ! -s "${STARTSTRUCTFILE}" ]; then
  die "${STARTSTRUCTFILE} does not exist or is empty!"
fi

shift 3

# process optional arguments
while [ -n "$1" ]; do
	case $1 in
		-rmin)
                  rmin=$2
                  shift
                  ;;
                -rmax)
                  rmax=$2
                  shift
                  ;;
                -numCPU)
                  numCPU=$2
                  shift
                  ;;
                -n)
                  n=$2
                  shift
                  ;;
                -outDir)
                  outDir=$2
                  shift
                  ;;
                -locusName)
                  locusName=$2
                  shift
                  ;;
                -?*)
                  echo "Unknown option (ignored): $1" >&2
                  ;;
        esac

        shift
done


# output options
OUTDIR=${outDir:-HiPR_output} # output folder
LOCUSNAME=${locusName:-UnnamedLocus} # locus name

# MCMC run parameters
NITER=${n:-100000}
NCPU=${numCPU:-16}
MODELTYPE=${modelType:-2} # set to 2
MINREADLEN=${rmin:-15}
MAXREADLEN=${rmax:-60}

LOGDIR="${OUTDIR}/logs"

# make output directories if do not exist already
mkdir -p "${OUTDIR}"
mkdir -p "${LOGDIR}"

NCHAIN=$(awk 'END{print NR/2}' ${STARTSTRUCTFILE}) # number of MCMC chains
head -n 1 ${STARTSTRUCTFILE} > ${OUTDIR}/${LOCUSNAME}.seq # primary sequence

echo "#chains=${NCHAIN}"

# set #CPUs to use
#ncpus=$(cat /proc/cpuinfo | awk '{if ($0~/^processor/) ++ncpu}END{print ncpu+0}') 
ncpus=$(get_ncpus)
echo >&2 "#cpus detected=${ncpus}"
echo >&2 "#cpus requested=${NCPU}"
if [ "${ncpus}" -gt "0" ] && [ "$NCPU" -gt "${ncpus}" ]; then
  echo >&2 "#requesting ${ncpus}"
  NCPU="${ncpus}"
fi

hasForks=$(perl -MForks::Super -e "print \"1\"" || echo "0")
[ "${hasForks}" -ne "1" ] && die "ERROR: Forks::Super module is not installed...\n"
#perl -MCPAN -e 'install Forks::Super'

perl ${HIPRMCMCwrapper} -id ${LOCUSNAME} -rmin ${MINREADLEN} -rmax ${MAXREADLEN} -burn_in 10 -samp 10 -t sd -numCPU ${NCPU} -n ${NITER} -m ${MODELTYPE} -outdir ${OUTDIR} ${COLLAPSEDREADFILE} ${MODRATEFILE} ${STARTSTRUCTFILE} > /dev/null && bash dms_mcmc_collect_structures.sh ${LOCUSNAME} ${OUTDIR} ${NCHAIN} > ${OUTDIR}/${LOCUSNAME}.mcmc.combined.struct && bash dms_get_majority_structure.sh ${OUTDIR}/${LOCUSNAME}.mcmc.combined.struct && cat ${OUTDIR}/${LOCUSNAME}.seq ${OUTDIR}/${LOCUSNAME}.mcmc.combined.struct.majority_consensus > ${OUTDIR}/${LOCUSNAME}.mcmc.majority_consensus


HiPR_structure=${OUTDIR}/${LOCUSNAME}.HiPR_structure.txt
HiPR_posterior=${OUTDIR}/${LOCUSNAME}.HiPR_posterior.txt

cp -p ${OUTDIR}/${LOCUSNAME}.mcmc.majority_consensus ${HiPR_structure} 
cp -p ${OUTDIR}/${LOCUSNAME}.mcmc.combined.struct.base_pairing_probs ${HiPR_posterior}

# concluding message
echo -e "\n\nHiPR algorithm complete\n------------------------------\n"
echo -e "Consensus secondary structure:\n${HiPR_structure}\n\n"
echo -e "Base-pairing posteriors at each nucleotide:\n${HiPR_posterior}\n"


: <<=cut
=pod

=head1 NAME

   HiPR.sh - High-throughput Probabilistic inference of RNA secondary structures. 

=head1 SYNOPSIS

   HiPR.sh reads_file rates_file structure_file [-outDir OutputDir] [-locusName LocusName] [-n Niter] [-numCPU Ncpus]  [-rmin MinReadLength] [-rmax MaxReadLength]"
   
   Mandatory arguments (need to be specified before optional arguments)
      reads_file -- File containing DMS-seq reads for RNA locus of interest (collapsed read format)
      rates_file -- File containing the initial estimates of per-nucleotide modification rates
      structure_file -- File containing the sequence and initial secondary structure of an RNA of interest

   Recognized optional command line arguments (need to be specified after mandatory arguments)
      -outDir <string>  -- Set name of output folder (default=HiPR_output)
      -locusName <string>  -- Set name of locus (default=UnnamedLocus)
      -n <integer> -- Set maximum number of MCMC iterations (default=100000)
      -numCPU <integer> -- Set number of CPUs to use (default=16)
      -rmin <integer> -- Set mininum read length (default=15)
      -rmax <integer> -- Set maximum read length (defauld=40)


=head1 DESCRIPTION

Estimate secondary structure and base pairing posteriors for a given RNA sequence based on the distribution of read fragments along the locus.

This program requires a file containing the sequence and initial secondary structure of an RNA of interest, a file containing DMS-seq reads, and a file containing the initial estimates of per-nucleotide modification rates. A Bayesian MCMC algorithm is then used to estimate the base pairing posterior that best fits the observed sequencing reads. The results and intermediate files are written to a directory (HiPR_output/ by default).


The output file HiPR_posterior.txt contains the base pairing posteriors at each nucleotide position, one entry per line.

The output file HiPR_structure.txt contains the consensus secondary structure.

=head2 Requirements

Mandatory files must be specified before any other optional arguments and must exist, otherwise a error or usage message will be shown. 

Software Requirements:
Linux-based operating system, or Mac OS.
Perl (v5.x). Important: Forks::Super module is required to run HiPR. To install, e.g., perl -MCPAN -e "install 'Forks::Super'", or see INSTALL.
Boost C++ libraries
Standard POSIX programs (awk, grep, bash)
C++ compiler (for re-compilation, if necessary)


=head1 LICENSE

MIT License https://opensource.org/licenses/MIT

Copyright (c) 2016 University of Pennsylvania

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

=head1 AUTHOR

Pavel Kuksa <pkuksa@upenn.edu>
Fan Li <fanli.gcb@gmail.com>
Li-San Wang <lswang@upenn.edu>

=cut
