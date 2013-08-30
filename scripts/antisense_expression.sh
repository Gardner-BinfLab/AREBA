
echo -e "Species\tEMBL\tExperiment_ID\tDataset_ID\tStrand_Correlation\tCore_CDS_On\tCore_Rfam_On\tAnnot_Concordance\tCoverage\tFraction_mapped_reads"
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
		#Coverage_mapped_reads=$(grep -E "^$species_directory.*$experiment_id.*$dataset_id.*${embl_name:0:8}" ../species_experiment_dataset_coverage_depth_mapped.stats | cut -f4,6);
		#if [ -z $Coverage_mapped_reads ]; then Coverage_mapped_reads=$(grep -E "^$species_directory.*$experiment_id.*$dataset_id" ../species_experiment_dataset_coverage_depth_mapped.stats | cut -f4,6); fi;

		#Correlation=$(echo `cat $plot | Rscript -e 'f=read.table("stdin"); cor(f[,1],f[,2]);'` | awk '{print $2}');
		if [ ! -f /home/suu13/projects/antisense/core_fasta/$embl_name*/$embl_name.$plot.HMM_Core.expression ]; then gff_plot2median.py -g /home/suu13/projects/antisense/core_fasta/$embl_name*/$embl_name.HMM_Core.gff -p $plot -a 1> /home/suu13/projects/antisense/core_fasta/$embl_name*/$embl_name.$plot.HMM_Core.expression 2> /dev/null; fi;
		if [ ! /home/suu13/projects/antisense/core_fasta/$embl_name*/$embl_name.$plot-rfam-filtered.expression ]; then gff_plot2median.py -g /home/suu13/projects/antisense/core_fasta/$embl_name*/$embl_name-rfam-filtered.gff -p $plot -a 1> /home/suu13/projects/antisense/core_fasta/$embl_name*/$embl_name.$plot-rfam-filtered.expression 2> /dev/null; fi;

		#Core_CDS=$(cat $plot.HMM_Core.expression | awk 'BEGIN{FS="$"; total_c=0; treshold_c=0;}{total_c++; if($7>1) treshold_c++;}END{print treshold_c/total_c}');
		#Core_Rfam=$(cat $plot.Rfam.expression | awk 'BEGIN{FS="$"; total_c=0; treshold_c=0;}{total_c++; if($7>1) treshold_c++;}END{print treshold_c/total_c}');
		#echo -e "$species\t$embl_name\t$experiment_id\t$dataset_id\t$Correlation\t$Core_CDS\t$Core_Rfam\t$Concordance\t$Coverage_mapped_reads";
	fi;
	done
cd /home/suu13/projects/areba/Quality_Score;

done
