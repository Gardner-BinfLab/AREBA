#iterate over directories and calculate cumulative stat files for ranking

for i in `ls */*embl | cut -d/ -f1`;

do
	
	cd $i;
	if [ -f max.plot ];
	then	
		genome_file=$(echo `ls *embl`| cut -d. -f1);
		echo $genome_file;
		#cat `ls *gff | grep -v "CDS\|gene\|annot"` | grep "RUF\|Rfam\|Pfam" | sponge $genome_file.cumulative.gff;
		#cat *PfamA*gff *rfam*gff *RUF*gff > $genome_file.cumulative.gff;
		cat ~/genomes/Pfam/*/*$genome_file* > ~/tmp/tmp.gff; #combine PfamAB	
		if [ -f max.plot ]; then gff_plot2median.py -g ~/tmp/tmp.gff -g ../../rfam_stinus_gff/$genome_file*Rfam.tblout.gff -g *RUF.gff -p max.plot > $genome_file.PfamAB.cumulative.stat; fi;
		cat $genome_file.PfamAB.cumulative.stat | sort -t$ -k4nr | sponge $genome_file.PfamAB.cumulative.stat;
		rm ~/tmp/tmp.gff;
	fi
	cd ..;
done
	
for i in `ls */*embl`;

do
	
	genome_file=$(echo $i | cut -d/ -f2 | cut -d. -f1);
	cd `echo $i | cut -d/ -f1`;
	if [ -f $genome_file-max.plot ];
	then	
		echo $genome_file;
		#cat `ls *gff | grep -v "CDS\|gene\|annot"` | grep "RUF\|Rfam\|Pfam" | sponge $genome_file.cumulative.gff;
		#cat *PfamA*gff *rfam*gff *RUF*gff > $genome_file.cumulative.gff;
		cat ~/genomes/Pfam/*/*$genome_file* > ~/tmp/tmp.gff;	
		if [ -f $genome_file-max.plot ]; then gff_plot2median.py -g ~/tmp/tmp.gff -g ../../rfam_stinus_gff/$genome_file*Rfam.tblout.gff -g $genome_file*RUF.gff -p $genome_file-max.plot > $genome_file.PfamAB.cumulative.stat; fi;
		cat $genome_file.PfamAB.cumulative.stat | sort -t$ -k4nr | sponge $genome_file.PfamAB.cumulative.stat;
		rm ~/tmp/tmp.gff;
	fi
	cd ..;
done


#table maker part
#for i in `cat Conservation_Counts_UPDATED.txt | cut -f1 | awk 'BEGIN{FS="_"}{print $1".*"$2".*"$3"|"$1".*"$3".*"$2}'`; do grep -E $i */*cumulative.stat; done


#for i in `cat Conservation_Counts_UPDATED.txt | cut -f1 | cut -d '_' --output-delimiter=".*" -f1,2,3`; do grep -E -n $i */*cumulative.stat; done
	
	
