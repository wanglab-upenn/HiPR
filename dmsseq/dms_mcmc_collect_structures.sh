LOCUSNAME=$1 
WDIR=$2 # directory with mcmc results
NCHAINS=${3:-100} # no chains to integrate
DOPROBPERITER=${4:-0} # set to 1 to output cumulative prob estimates

basefn=${WDIR}/${LOCUSNAME}.mcmc # prefix for MCMC results
for i in `seq 1 ${NCHAINS}`; do
  #echo "${i}"
  chain=$(printf "%03d" ${i}) # chain number in xxx format
  locus_mcmc_file=${basefn}.${chain}.txt
  locus_log_file=${basefn}.${chain}.log
  locus_sampled_struct_file=${basefn}.${chain}.mcmc_struct
  locus_prob_per_iter_file=${basefn}.${chain}.prob_per_iter
  locus_accepted_samples_file=${basefn}.${chain}.accepted_struct
  paste <(tail -n +2 ${locus_mcmc_file}) ${locus_log_file} | \
   awk 'BEGIN{ ll=-100000;
               sampFile="'${locus_sampled_struct_file}'";
               probFile="'${locus_prob_per_iter_file}'";
               acceptedFile="'${locus_accepted_samples_file}'";
               doProbEstimates="'${DOPROBPERITER}'";
               OFS="\t";
               chain=("chain_" "'${chain}'")
             }
        { 
       iter=$7;
       if ($8 !~ /accepted/) next; # get only "accepted" structures
       # output accepted structures
       #print $1, $2 > acceptedFile;
       ll_current=$2;
       if (ll_current>ll)
       {
         nsamp+=1;
         ll=ll_current;
         structure=$1;
         if (doProbEstimates==1)
         {
           n=split(structure,s,"");
           for (j=1;j<=n;j++)
           {
             a=s[j] #substr(structure,j,1);
             if (a=="(" || a==")") probs[j]++;
           }
           #print current base-pairing prob estimates
           prob_current=(probs[1]+0.0)/nsamp
           for (j=2;j<=n;j++)
             prob_current=(prob_current "\t" (probs[j]+0.0)/nsamp)
           printf "%s\n", prob_current > probFile
         }
         print structure, ll_current, chain # > sampFile
       }
     }'
done | \
 sort -k2,2nr 


# combine samples from different chains
# samples will be sorted in the decreasing order of LL
#combinedStructFile=${basefn}.combined.mcmc_struct
#>${combinedStructFile}
#cat ${basefn}.*.mcmc_struct | sort -k2,2nr > ${combinedStructFile}


