
echo -e "#RfamID\tRfamAlias\tEMBLAccession\tExpressionRank"



for file in `ls */*$1`;
do
embl_accession=$(echo $file | cut -d/ -f2 | cut -d. -f1);
cat $file | gawk 'match($0,/ID=(.*);Alias=(.*);E-value/,alias)  {print alias[1],"\t",alias[2],"\t","'$embl_accession'","\t",NR}';

done




: << 'END'
for keyword in `cat */*rfam*.gff | gawk 'match($0,/Alias=(.*);/,m) {print m[1]}' | sort -n | uniq`;
do

	order=$(for i in `grep -E Alias=$keyword */*cumulative.stat | cut -f1 | cut -d: -f2 | cut -d. -f1`; do grep -n $i ../order.embl_2.txt; done | sort -n | uniq | wc -l); #uniq ile cozdum redundancy meselesini
	grep -E -n -m1 Alias=$keyword */*cumulative.stat | gawk 'BEGIN{FS=":"}{match($1,/\/(.*)\.cumulative/,m); printf"Rfam_%s_%s\t%s\t\t\t\t\t%s\n","'$keyword'",m[1],$2,"'$order'";}'

	
done

END
