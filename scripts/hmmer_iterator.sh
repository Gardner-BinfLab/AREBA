

for i in `ls ../Quality_Score/*embl | cut -d/ -f3 | awk '{print substr($0,1,length($0)-5)}'`; #iterate over embl names

do
	
	if [ -f ../Quality_Score/$i-CDS-AA.fasta ] && [ ! -f ../Quality_Score/$i.HMM_Core.gff ]; #check CDS fasta file available and hmm gff available
	then
		echo $i;
		hmmsearch -o /dev/null -E 0.05 --tblout ../Quality_Score/$i.tblout.txt --domtblout ../Quality_Score/$i.domtblout.txt /home/suu13/projects/areba/core_bacteria_archea/combined.bact.arch.hmm ../Quality_Score/$i-CDS-AA.fasta; #do an hmm search based on combined profiles of bacteria and archea on all CDS fasta files of genomes
		grep -v ^# ../Quality_Score/$i.tblout.txt | awk 'BEGIN{ad=foo}{if($3!=ad) {ad=$3; print $0}}' > ../Quality_Score/$i.tblout.tmp; #select best scoring uniq matches
		hmm_Parser.py -t ../Quality_Score/$i.tblout.tmp > ../Quality_Score/$i.HMM_Core.gff; #parse tmp and convert it to gff
		rm ../Quality_Score/$i.tblout.tmp; #remove tmp file
	fi
done
