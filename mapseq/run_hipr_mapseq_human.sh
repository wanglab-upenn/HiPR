source "dmsmapseq.ini"

INBAM=$1 # prepared BAM file
LOCIBED=$2 # list of loci, 1-based coord!!! BED6 format
OUTDIR=$3 # output directory
collectOnly=${4:-0} # 1=only run post-analysis (collection, consensus, etc)

GENOMEFA=/project/wang/pkuksa/datasets/hg19/hg19.fa
INITRATES=`dirname $0`/gold.rates.dmsmapseq.ACUG

numChains=100 #100
maxIter=10000 #100000
numCPU=16
modelType=2
minReadLength=45
maxReadLength=60
WITHINONLY=1
MAX_MISM_CNT=5
OPTTYPE="sd" # s=struct only; sd=struct+params

BINDIR=/project/wang/pkuksa/bin

export RNASEQFOLD="perl -I /appl/perl-5.18.2/lib/site_perl/5.18.2/ ${BINDIR}/RNASeqFold/HiPR_dmsmapseq_prod.pl"

echo "bash ${SCRIPTDIR}/run_dmsmapseq_pipeline.sh ${INBAM} ${GENOMEFA} ${LOCIBED} ${INITRATES} ${OUTDIR} ${WITHINONLY} ${numChains} ${maxIter} ${numCPU} ${modelType} ${minReadLength} ${maxReadLength} ${collectOnly} ${MAX_MISM_CNT} ${OPTTYPE}"
bash ${SCRIPTDIR}/run_dmsmapseq_pipeline.sh ${INBAM} ${GENOMEFA} ${LOCIBED} ${INITRATES} ${OUTDIR} ${WITHINONLY} ${numChains} ${maxIter} ${numCPU} ${modelType} ${minReadLength} ${maxReadLength} ${collectOnly} ${MAX_MISM_CNT} ${OPTTYPE}


