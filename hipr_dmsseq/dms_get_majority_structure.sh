INSTRUCT=$1 # input structures, one structure per line
numSamples=$(wc -l ${INSTRUCT} | awk '{print $1}')
structLength=$(awk '{print length($1); exit}' ${INSTRUCT});
awk 'function top(a)        {return a[a[0]]}
     function push(a,x,i) {i=++a[0]; a[i]=x; return i}
     function pop(a,x,i)
     {
        i=a[0]--;  
        if (!i) {return -1} else {x=a[i]; delete a[i]; return x}
     }
     BEGIN{ sep = "\t"; OFS=sep; }
     {
       len=split($1,str,"");
       split("",stack,"");
       for (j=1;j<=len;++j)
       {
         a=str[j];
         if (a=="(") push(stack,j,i);
         if (a==")")
         {
           bp=pop(stack); 
           #print "pair:",bp,j;
           bpC[(bp sep j)]++;
         }
       }
     }END{for (bp in bpC) print bp, bpC[bp];}' ${INSTRUCT} | \
       sort -k3,3nr  > ${INSTRUCT}.bp_freq #| \

     awk 'BEGIN{OFS="\t"; m='${numSamples}'+0;}
          { 
            bpcount = $3;
            bpprob = bpcount / m
            print $1, $2, $3, bpprob
          }' ${INSTRUCT}.bp_freq > ${INSTRUCT}.bp_prob

     awk 'BEGIN{OFS="\t"; m='${numSamples}'+0; n='${structLength}'+0;}
          { 
            bpcount = $3;
            bpprob = bpcount / m;
            bp_i = $1;
            bp_j = $2;
            b[bp_i]+=bpcount
            b[bp_j]+=bpcount
            #print $1, $2, $3, bpprob
          }
          END{ for (i=1; i<=n; ++i)
                #printf "%d\t%.6f\t%d\t%d\n", b[i], b[i] / m, m, n
                printf "%.6f\n", b[i] / m
             }' ${INSTRUCT}.bp_freq > ${INSTRUCT}.base_pairing_probs

     
     awk 'BEGIN{m='${numSamples}'+0;maj=int(m/2.0+1);}
          {
            # base-pairs
            bpcount = $3;
            if (bpcount >= maj)
              print;
          }' ${INSTRUCT}.bp_freq > ${INSTRUCT}.majority_bp

     awk 'BEGIN{len='${structLength}'+0; for (i=1;i<=len;++i) scons[i]=".";}
          {
            scons[$1+0]="(";
            scons[$2+0]=")";
             
          }END{for (i=1;i<=len;i++) printf "%s", scons[i]; printf "%s","\n";}'  ${INSTRUCT}.majority_bp > ${INSTRUCT}.majority_consensus
     awk '{len=split($1,s,"");
           # binary base-paring state vector
           for (i=1;i<=len;++i)
           {
             a = s[i];
             if (a==".") printf "%s\t", "0";
             else printf "%s\t", "1";
           }
          }' ${INSTRUCT}.majority_consensus > ${INSTRUCT}.majority_binary
