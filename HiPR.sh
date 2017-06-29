# functions
function die {
   echo "$1" >&2
   exit 1
}

# check HiPR executables
scriptDir=`dirname $0`
PATH=$PATH:${scriptDir}
mcmcCmd=HiPR_MCMC
mcmcCmdWrapper=HiPR_MCMC_wrapper.pl
HIPRMCMC=${HIPRMCMC:-$(type -P ${mcmcCmd})}
command -v "${HIPRMCMC}" > /dev/null 2>&1 || { echo >&2 "$0 requires ${mcmcCmd}, but it is not found. Please set HIPRMCMC or make sure ${mcmcCmd} is in the path"; exit 1; }
HIPRMCMCwrapper="${scriptDir}/${mcmcCmdWrapper}"
if [ ! -s "${HIPRMCMCwrapper}" ]; then 
  echo >&2 "$0 requires ${mcmcCmdWrapper}, but it is not found. Please make sure ${mcmcCmdWrapper} is in the path";
  exit 1;
fi

# required inputs [must be first 3 arguments]
COLLAPSEDREADFILE=$1 # DMS.collapsed.reads
MODRATEFILE=$2 # per-nucleotide modification rates for unpaired/paired bases
STARTSTRUCTFILE=$3 # e.g., structures sampled from Boltzmann distribution

if [ $# -lt 3 ]; then
  #echo >&2 "USAGE: $0 reads_file rates_file structure_file [OutputDir] [LocusName] [Niter] [Ncpus] [ModelType] [MinReadLength] [MaxReadLength]"
  echo >&2 "USAGE: $0 reads_file rates_file structure_file [-outDir OutputDir] [-locusName LocusName] [-n Niter] [-numCPU Ncpus]  [-rmin MinReadLength] [-rmax MaxReadLength]"
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
ncpus=$(cat /proc/cpuinfo | awk '{if ($0~/^processor/) ++ncpu}END{print ncpu+0}') 
echo >&2 "#cpus detected=${ncpus}"
echo >&2 "#cpus requested=${NCPU}"
if [ "${ncpus}" -gt "0" ] && [ "$NCPU" -gt "${ncpus}" ]; then
  echo >&2 "#requesting ${ncpus}"
  NCPU="${ncpus}"
fi

perl ${HIPRMCMCwrapper} -id ${LOCUSNAME} -rmin ${MINREADLEN} -rmax ${MAXREADLEN} -burn_in 10 -samp 10 -t sd -numCPU ${NCPU} -n ${NITER} -m ${MODELTYPE} -outdir ${OUTDIR} ${COLLAPSEDREADFILE} ${MODRATEFILE} ${STARTSTRUCTFILE} > /dev/null && bash dms_mcmc_collect_structures.sh ${LOCUSNAME} ${OUTDIR} ${NCHAIN} > ${OUTDIR}/${LOCUSNAME}.mcmc.combined.struct && bash dms_get_majority_structure.sh ${OUTDIR}/${LOCUSNAME}.mcmc.combined.struct && cat ${OUTDIR}/${LOCUSNAME}.seq ${OUTDIR}/${LOCUSNAME}.mcmc.combined.struct.majority_consensus > ${OUTDIR}/${LOCUSNAME}.mcmc.majority_consensus
