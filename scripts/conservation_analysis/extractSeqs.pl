use strict;
use warnings;
use Bio::SeqIO;
use Data::Dumper;
use Cwd;
use File::Find;
my $dir = cwd;

#SYNTAX: perl extractSeqs.pl > /embl.new_genomes/extracted_seqs.fa
#This code reads the coordinates in rufEntries_clean.gff and extract sequences from relevant genome directories in fastA format.
#biohisham@gmail.com

open (my $coordFH, "<","$dir/embl.new_genomes/rufEntries_clean.gff") or die ("could not open file $!\n");
while(my $line = <$coordFH>){
	chomp $line;
	#print $line;
	my ($id, $coord1, $coord2) = split("\t",$line);
	if($coord1 < $coord2){	#forward strand
		extractForward($id, $coord1, $coord2);	
		}elsif($coord1 > $coord2){	#reverse strand
			extractBackward($id, $coord1, $coord2);		
			}
	}


sub extractForward{
	my $id = shift;
	my $coord1 = shift;
	my $coord2 = shift;
	my $fileName = getFileName($id);
	print STDERR $fileName,$/;
	my $in = Bio::SeqIO->new(-file=>$fileName, -format=>'embl');
	while(my $seq = $in->next_seq){
		print  ">","$id","_","$coord1","_", "$coord2","\n";
		print $seq->subseq($coord1, $coord2), "\n";
		}
	}

sub extractBackward{
	my $id = shift;
	my $coord1 = shift;
	my $coord2 = shift;
	my $fileName = getFileName($id);
	print STDERR $fileName, $/;
	my $in = Bio::SeqIO->new(-file=>$fileName, -format=>'embl');
	while(my $seq = $in->next_seq){
		print ">$id","_","$coord2","_", "$coord1","\n";
		print  $seq->subseq($coord2, $coord1), "\n";
		}
	}

sub getFileName{
	our $genomeID = shift;
	our $file; #variable to hold File::Find object
	find(\&wanted,"$dir/embl.new_genomes");
	sub wanted{
		if($File::Find::name =~ /$genomeID.*?\.embl$/){
			$file = $File::Find::name;			
			}
		}
	return $file;	
	}

##cat *.embl>AllGenomes.embl
#make a directory called nhmmer
#mv AllGenomes.embl nhmmer
#mv embl.new_genomes/extracted_seqs.fa embl.new_genomes/nhmmer
#run runNHMMER.pl
