

query=





"5S_rRNA\|tmRNA\|SSU_rRNA\|Bacteria_small_SRP\|RNaseP\|Bacteria_large_SRP"

grep $RNAon $qstatrfam | sort -t$ -k8nr | cut -d'$' -f8 | head -n 10 | awk 'BEGIN{i=0; result=0;}{array[i]=$1; i++;}END{for(a=0;a<i;a++) if(array[a]>=0.5) result+=1; printf"ncRNA On : %f\n",result/i;}'


grep "ribosomal" AE017225.1.gff.qstat | sort -t$ -k8nr | cut -d'$' -f8 | head -n 10 | awk 'BEGIN{i=0; result=0;}{array[i]=$1; i++;}END{for(a=0;a<i;a++) if(array[a]>=0.5) result+=1; print result/i;}'


grep "dnaG\|rplA\|rpmA\|rpsS\|pyrG\|rpsM\|rpsJ\|rpsM\|rpsC\|rplC\|rplB\|rpsB\|rpsE" AE017225.1.gff.qstat | sort -t$ -k8nr | cut -d'$' -f8 | head -n 10 | awk 'BEGIN{i=0; result=0;}{array[i]=$1; i++;}END{for(a=0;a<i;a++) if(array[a]>=0.5) result+=1; print result/i;}'


grep CDS GSE13543-max.gff.qstat | sort -t$ -k6nr | cut -d'$' -f2,6 --output-delimiter=" " | head -n 100 | awk '{if($2/$1>2) result+=1}END{printf"Concordance:%f\n",result/100;}'

cat max.plot | Rscript -e 'f=read.table("stdin"); cor(f[,1],f[,2]);'
