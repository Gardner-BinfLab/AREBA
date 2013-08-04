
echo -e "EMBL\tExperiment_Name\tCorrelation\tCore_On\tCore_Rfam_On\tConcordance"
cat embl.directory.map.txt | while read i;
do

embl_name=$(echo $i | cut -d' ' -f2;);
species_directory=$(echo $i | cut -d' ' -f1;);
cd plot/$species_directory;


	for plot in `ls */*stranded.plot`;
	do
	name_of_experiment=$(echo $plot | cut -d. -f1);
	correlation=$(echo `cat $plot | Rscript -e 'f=read.table("stdin"); cor(f[,1],f[,2]);'` | awk '{print $2}');
	gff_plot2median.py -g /home/suu13/projects/areba/Quality_Score/$embl_name.HMM_Core.gff -p $plot -a 1> $plot.CDS.expression 2> /dev/null;
	gff_plot2median.py -g /home/suu13/projects/areba/Quality_Score/$embl_name.rfam.core.gff -p $plot -a 1> $plot.Rfam.expression 2> /dev/null;
	gff_plot2median.py -g /home/suu13/projects/areba/Quality_Score/$embl_name.all.gff -p $plot 1> $plot.all.expression 2> /dev/null;
	Core_CDS=$(cat $plot.CDS.expression | awk 'BEGIN{FS="$"; total_c=0; treshold_c=0;}{total_c++; if($7>1) treshold_c++;}END{print treshold_c/total_c}');
	Core_Rfam=$(cat $plot.Rfam.expression | awk 'BEGIN{FS="$"; total_c=0; treshold_c=0;}{total_c++; if($7>1) treshold_c++;}END{print treshold_c/total_c}');
	Concordance=$(cat $plot.all.expression | awk 'BEGIN{FS="$"; total_c=0; treshold_c=0;}{total_c++; if($4>$6) treshold_c++;}END{print treshold_c/total_c}');
	echo -e "$embl_name\t$name_of_experiment\t$correlation\t$Core_CDS\t$Core_Rfam\t$Concordance";
	done
cd /home/suu13/projects/areba/Quality_Score;

done
