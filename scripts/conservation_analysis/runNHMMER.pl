use strict;
use warnings;
use Cwd;
use Data::Dumper;
use Bio::SeqIO;

#This program extracts each sequence from extracted_seqs.fa into a file named after the FastA ID for that sequence.
#The extracted sequence is searched against AllGenomes.embl and the results are saved into nhmmer/$FastAID
#SYNTAX: Perl runNHMMER.pl
#biohisham@gmail.com
my @fileName_array; #at the end move over the array to concatenate output files based on seqName
my $dir = getcwd;
#print $dir,$/;
my $in = Bio::SeqIO->new(-file=> "$dir/embl.new_genomes/nhmmer/extracted_seqs.fa", -format=>'fasta');
while(my $seq= $in->next_seq){
	my $fileName= $seq->id; #the FastA_ID based fileName in directory nhmmer.
	print "creating FastA file\n";
	print $fileName, "\n";
	push @fileName_array, $fileName; 
	my $out = Bio::SeqIO->new(-file=>">$dir/embl.new_genomes/nhmmer/$fileName", -format=>'fasta');
	$out->write_seq($seq);
	#call nhmmer on the file
	#nhmmer -o nhmmer_output --tblout nhmmer_results_table -E 0.001 --rna --acc --toponly testgeneSeq.fa AllGenomes.embl
	print "Running nhmmer \n";
	#nhmmer arguments definition
	my $outFile = "$dir/embl.new_genomes/nhmmer/$fileName"."_"."out";
	my $tabulatedoutput="$dir/embl.new_genomes/nhmmer/$fileName"."_"."table";
	next if (-e $outFile && -e $tabulatedoutput);
	my $eval = "0.001";
	my $query = "$dir/embl.new_genomes/nhmmer/$fileName";
	my $database = "$dir/embl.new_genomes/nhmmer/AllGenomes.embl";
	system("nhmmer -o $outFile --tblout $tabulatedoutput -E $eval --rna --acc --toponly $query $database");
	}

