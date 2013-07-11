#iterate over directories and calculate cumulative stat files for ranking


for i in `ls */*embl | cut -d/ -f1`;

do
	cd $i;
	if [ -f max.plot ];
	then	
		genome_file=$(echo `ls *embl`| cut -d. -f1);
		cat `ls *gff | grep -v "gene\|annot"` | grep "CDS\|RUF\|Rfam" | sponge $genome_file.cumulative.gff;
		gff_plot2median.py -g $genome_file.cumulative.gff -p max.plot > $genome_file.cumulative.stat;
		cat $genome_file.cumulative.stat | sort -t$ -k4nr | sponge $genome_file.cumulative.stat;
	fi
	cd ..;
done
	
	
	
