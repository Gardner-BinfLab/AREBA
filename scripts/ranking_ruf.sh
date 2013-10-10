

echo -e "#RUF_ID\tStrand\tEMBLAccession\tExpressionRank"

for file in `ls */*$1`; #use extensin of file as an argument
do

cat $file | gawk 'match($0,/misc_RNA.*RUF/,alias)  {print "RUF_"$4"_"$5,"\t",$7,"\t",$1,"\t",NR}';

done
