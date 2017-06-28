
# check HiPR
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

# required inputs
COLLAPSEDREADFILE=$1 # DMS.collapsed.reads
STARTSTRUCTFILE=$2 # e.g., structures sampled from Boltzmann distribution
MODRATEFILE=$3 # per-nucleotide modification rates for unpaired/paired bases

if [ $# -lt 3 ]; then
  echo >&2 "USAGE: $0 reads_file rates_file structure_file"
  exit 1;
fi

# output options
OUTDIR=${4:-HiPR_output} # output folder
LOCUSNAME=${5:-UnnamedLocus} # locus name

# MCMC run parameters
NITER=${6:-100000}
NCPU=${7:-16}
MODELTYPE=${8:-2} # set to 2
MINREADLEN=${9:-15}
MAXREADLEN=${10:-60}


LOGDIR="${OUTDIR}/logs"


#QSUB="bsub"
#jobName="mcmc_${LOCUSNAME}"

mkdir -p ${OUTDIR} # make output dir if does not exist already
mkdir -p "${LOGDIR}"

NCHAIN=$(awk 'END{print NR/2}' ${STARTSTRUCTFILE}) # number of MCMC chains
head -n 1 ${STARTSTRUCTFILE} > ${OUTDIR}/${LOCUSNAME}.seq # primary sequence

echo "#chains=${NCHAIN}"

ncpus=$(cat /proc/cpuinfo | awk '{if ($0~/^processor/) ++ncpu}END{print ncpu+0}') 
  echo >&2 "#cpus detected=${ncpus}"
  echo >&2 "#cpus requested=${NCPU}"
if [ "${ncpus}" -gt "0" ] && [ "$NCPU" -gt "${ncpus}" ]; then
  echo >&2 "#requesting ${ncpus}"
  NCPU="${ncpus}"
fi
perl ${HIPRMCMCwrapper} -id ${LOCUSNAME} -rmin ${MINREADLEN} -rmax ${MAXREADLEN} -burn_in 10 -samp 10 -t sd -numCPU ${NCPU} -n ${NITER} -m ${MODELTYPE} -outdir ${OUTDIR} ${COLLAPSEDREADFILE} ${MODRATEFILE} ${STARTSTRUCTFILE} > /dev/null && bash dms_mcmc_collect_structures.sh ${LOCUSNAME} ${OUTDIR} ${NCHAIN} > ${OUTDIR}/${LOCUSNAME}.mcmc.combined.struct && bash dms_get_majority_structure.sh ${OUTDIR}/${LOCUSNAME}.mcmc.combined.struct && cat ${OUTDIR}/${LOCUSNAME}.seq ${OUTDIR}/${LOCUSNAME}.mcmc.combined.struct.majority_consensus > ${OUTDIR}/${LOCUSNAME}.mcmc.majority_consensus
