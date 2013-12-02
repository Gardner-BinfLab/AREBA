
echo -e "Species\tEMBL\tExperiment_ID\tDataset_ID\tStrand_Correlation\tCore_CDS_On\tCore_Rfam_On\tPPV\tSpecificity\tConcordance_New_Metric\tCoverage\tFraction_mapped_reads"
cat embl.directory.map.txt | while read i;
do

embl_name=$(echo $i | cut -d' ' -f2;);
species_directory=$(echo $i | cut -d' ' -f1;);
cd plot/$species_directory;


	for plot in `ls */*stranded.plot`;
	do
	if [ `head /home/suu13/projects/areba/Quality_Score/$embl_name.embl | gawk 'match($0,/; ([0-9]*) /,m) {print m[1]}'` -eq `wc -l $plot | cut -d' ' -f1` ]; #embl dosyasi ve plot dosyasinin uyumlu oldugunu kontrol etmek icin.
	then
		experiment_id=$(echo $plot | cut -d/ -f1);
		dataset_id=$(echo $plot | cut -d. -f1 | cut -d/ -f2);
		species=$(grep -m1 $dataset_id ../species_experiment_dataset_coverage_depth_mapped.stats | cut -f1);
		Coverage_mapped_reads=$(grep -E "^$species_directory.*$experiment_id.*$dataset_id.*${embl_name:0:8}" ../species_experiment_dataset_coverage_depth_mapped.stats | cut -f4,6);
		if [ -z $Coverage_mapped_reads ]; then Coverage_mapped_reads=$(grep -E "^$species_directory.*$experiment_id.*$dataset_id" ../species_experiment_dataset_coverage_depth_mapped.stats | cut -f4,6); fi;

		Correlation=$(echo `cat $plot | Rscript -e 'f=read.table("stdin"); cor(f[,1],f[,2]);'` | awk '{print $2}');
		if [ ! -f $plot.CDS.expression ]; then gff_plot2median.py -g /home/suu13/projects/areba/Quality_Score/$embl_name.HMM_Core.gff -p $plot -a 1> $plot.CDS.expression 2> /dev/null; fi;
		if [ ! -f $plot.Rfam.expression ]; then gff_plot2median.py -g /home/suu13/projects/areba/Quality_Score/$embl_name.rfam.core.gff -p $plot -a 1> $plot.Rfam.expression 2> /dev/null; fi;
		if [ ! -f $plot.all.expression ]; then gff_plot2median.py -g /home/suu13/projects/areba/Quality_Score/$embl_name.all.gff -p $plot 1> $plot.all.expression 2> /dev/null; fi;
		#Genome_Median=$(gff_plot2median.py -g foo -p $plot -m | cut -d: -f2); #calculate genomic median
		#if [ $Genome_Median -eq 0 ]; then let $Genome_Median=1; fi;  #if genomic median is 0, make it 1
		Genome_Median=10;
		Core_CDS=$(cat $plot.CDS.expression | awk 'BEGIN{FS="$"; total_c=0; treshold_c=0;}{total_c++; if($4>='"$Genome_Median"') treshold_c++;}END{print treshold_c/total_c}'); #check core CDS on
		Core_Rfam=$(cat $plot.Rfam.expression | awk 'BEGIN{FS="$"; total_c=0; treshold_c=0;}{total_c++; if($4>='"$Genome_Median"') treshold_c++;}END{print treshold_c/total_c}'); #check selected rfam on
		#Concordance=$(Concordance.py -p $plot -g /home/suu13/projects/areba/Quality_Score/$embl_name.gff3 -g /home/suu13/projects/areba/Quality_Score/$embl_name.all.gff | awk '{print $2}'); #check concordance of whole genome for all annotations
		PPV_Spec=$(Concordance.py -a -p $plot -g /home/suu13/projects/areba/Quality_Score/$embl_name.gff3 -g /home/suu13/projects/areba/Quality_Score/$embl_name.all.gff | cut -d: -f2); #do not forget to put -a
		#echo -e "$species\t$embl_name\t$experiment_id\t$dataset_id\t$Correlation\t$Core_CDS\t$Core_Rfam\t$Concordance\t$Coverage_mapped_reads";
		echo -e "$species\t$embl_name\t$experiment_id\t$dataset_id\t$Correlation\t$Core_CDS\t$Core_Rfam\t$PPV_Spec\t$Coverage_mapped_reads";

	fi;
	done
cd /home/suu13/projects/areba/Quality_Score;

done
