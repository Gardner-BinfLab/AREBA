
echo -e "#hit\tExpression_Rank\t\t\t\t\torders_seen"


for keyword in `cat */*rfam*.gff | gawk 'match($0,/Alias=(.*);/,m) {print m[1]}' | sort -n | uniq`;
do

	order=$(for i in `grep -E Alias=$keyword */*cumulative.stat | cut -f1 | cut -d: -f2 | cut -d. -f1`; do grep -n $i ../order.embl_2.txt; done | sort -n | uniq | wc -l); #uniq ile cozdum redundancy meselesini
	grep -E -n -m1 Alias=$keyword */*cumulative.stat | gawk 'BEGIN{FS=":"}{match($1,/\/(.*)\.cumulative/,m); printf"Rfam_%s_%s\t%s\t\t\t\t\t%s\n","'$keyword'",m[1],$2,"'$order'";}'

	
done
