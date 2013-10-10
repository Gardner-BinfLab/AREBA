echo -e "#CDS_ID\tStrand\tEMBLAccession\tExpressionRank"

for file in `ls */*$1`;
do
cat $file | gawk 'match($0,/EMBL.*CDS/,alias)  {print "CDS_"$4"_"$5,"\t",$7,"\t",$1,"\t",NR}';

done
