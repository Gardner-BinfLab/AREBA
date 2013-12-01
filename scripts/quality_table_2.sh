echo -e "EMBL\tExperiment_Name\tCorrelation\tCore_On\tCore_Rfam_On\tConcordance\tCoverage"
cat plot/species_experiment_dataset_coverage_depth_mapped.stats | while read i;
do
	
	experiment_id=$(echo $i | awk '{print $2}');
	data_set_id=$(echo $i | awk '{print $3}');
	species_name=$(echo $i | awk '{print $1}');
	coverage=$(echo $i | awk '{print $4}');
	fraction=$(echo $i | awk '{print $6}');


		cd plot/$species_name;
			#echo $data_set_id;
		if [ -f */$data_set_id.stranded.plot ];
			then	
			plot=$(echo `ls */$data_set_id.stranded.plot`);
	
			Correlation=$(echo `cat $plot | Rscript -e 'f=read.table("stdin"); cor(f[,1],f[,2]);'` | awk '{print $2}');
			if [ ! -f $plot.CDS.expression ]; then gff_plot2median.py -g /home/suu13/projects/areba/Quality_Score/$embl_name.HMM_Core.gff -p $plot -a 1> $plot.CDS.expression 2> /dev/null; fi;
			#if [ ! -f $plot.Rfam.expression ]; then gff_plot2median.py -g /home/suu13/projects/areba/Quality_Score/$embl_name.rfam.core.gff -p $plot -a 1> $plot.Rfam.expression 2> /dev/null; fi;
			#if [ ! -f $plot.all.expression ]; then gff_plot2median.py -g /home/suu13/projects/areba/Quality_Score/$embl_name.all.gff -p $plot 1> $plot.all.expression 2> /dev/null; fi;
			#Core_CDS=$(cat $plot.CDS.expression | awk 'BEGIN{FS="$"; total_c=0; treshold_c=0;}{total_c++; if($7>1) treshold_c++;}END{print treshold_c/total_c}');
			#Core_Rfam=$(cat $plot.Rfam.expression | awk 'BEGIN{FS="$"; total_c=0; treshold_c=0;}{total_c++; if($7>1) treshold_c++;}END{print treshold_c/total_c}');
			#Concordance=$(cat $plot.all.expression | awk 'BEGIN{FS="$"; total_c=0; treshold_c=0;}{total_c++; if($4>$6) treshold_c++;}END{print treshold_c/total_c}');
			echo -e "$species_name\t$experiment_id\t$data_set_id\t$Correlation\t$Core_CDS\t$Core_Rfam\t$Concordance\t$coverage";
		fi;
		cd /home/suu13/projects/areba/Quality_Score;

done
