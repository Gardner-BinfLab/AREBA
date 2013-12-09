

for i in `ls ../Quality_Score/*embl | cut -d/ -f3 | awk '{print substr($0,1,length($0)-5)}'`; #iterate over embl names

	do
	echo $i;
	embl2CDS_Translated.py -e ../Quality_Score/$i.embl > ../Quality_Score/$i-CDS-AA.fasta;

done


#hmmsearch -o /dev/null --tblout tblout.txt --domtblout domtblout.txt combined.hmm CP001472.1-CDS-AA.fasta
#cat tblout.txt | awk 'BEGIN{ad=foo}{if($3!=ad) {ad=$3; print $3,$5;}}'
