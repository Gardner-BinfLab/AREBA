use strict;
use warnings;
use Cwd;
use File::Find;
use Bio::SeqIO;
use Scalar::Util qw(looks_like_number);
use Data::Dumper;
my $dir = getcwd;

#print $dir,$/;
#Syntax "perl nhmmerTableParse_v1.pl > Conservation_Counts.txt"
#This program identifies if a nhmmer table hit was seen in more than one species.
#It collates information from different files generated upstream.
#It retains the hit ID, species counts, their IDs.
#It summarize these hits across different taxonomies using the information in "bacteria.taxonomy.ClassMatch.txt" and "bacteria.taxonomy.NonClass.txt".
#it returns the sequence queried using the relevant file to the nhmmer table output.
#biohisham@gmail.com
my %hash;
my ($seqFile,$tableFile);
print "#hit\tgenomes\ttotal_copy_number\tclasses_seen\tgenomes_in_class\tseq\n";
find(\&wanted, "$dir/embl.new_genomes/nhmmer");

sub wanted{
	if($File::Find::name =~ /(.*?)_table/){
		#print $File::Find::name,$/;
		$tableFile = $File::Find::name;
		$seqFile = $1;
		#print $seqFile, "\n";
		#print $tableFile, $/;;
		nhmmerTableSeqRead($seqFile,$tableFile);
		}
	}

sub nhmmerTableSeqRead{
	$seqFile = shift;
	$tableFile = shift;
#	print $tableFile,$/;
	open (my $fh,"<", "$tableFile") or die ("could not open file of tabulated nhmmer output $! \n");

	while(my $line = <$fh>){
		next if $line =~/^#/;
		#Skipping the organism name here cuz it has spaces that will interfere with 'spil'ling.
		#The other nhhmer table output columns may come in handy later
		my ($target, $acc1, $query, $acc2, $hmmFrom, $hmmTo, $aliFrom, $aliTo, $envFrom, $envTo, $seqLen, $strand, $e, $score, $bias) = split(/\s+/,$line);
		#for each $query count $target
		$hash{$query}{$target}++;	
		my $in = Bio::SeqIO->new(-file=>"$seqFile", -format =>'fasta');
		for (my $seq=$in->next_seq){
			 $hash{$query}{'seq'}=$seq->seq;
			}
		}


	}

	#Print query, number of hits, breakdown of hits per ID, seq.
	#See if hit is in bacteria.taxonomy.ClassMatch.txt or bacteria.taxonomy.NonClass.txt
	foreach my $key (keys %hash){
		print $key,"\t";
		my $genomes = keys (%{$hash{$key}}) - 1; #Subtract the 'seq' key from the others.
		print $genomes, "\t";
		#get the total copy numbers and the genome IDs to scan "bacteria.taxonomy.ClassMatch.txt" and "bacteria.taxonomy.NonClass.txt".
		my $totalCopies;
		while (my ($target,$value) = each (%{$hash{$key}})){
			#print checkClassed($target),"\t";
			#the return of subroutines checkClassed determines seen or not.
			$hash{$key}{'classified'}{checkClassed($target)}++;
			looks_like_number($value) ? $totalCopies +=$value : next;
			}
		print $totalCopies,"\t";
		my $taxaCount=0; #number of taxa seen for the keys.
		my $genomeCount=0; #number of genomes seen in each tax group
		for my $tax (keys %{$hash{$key}{'classified'}}){
			#print $tax,$/;
			#number of classes
			next if $tax eq "";
			++$taxaCount;
			$genomeCount += $hash{$key}->{'classified'}->{$tax};
			}
		#return classes seen (total of times seen in each class is equal to number of genomes)
		print $taxaCount,"\t"; #Classes seen
		#total genomes in class
		print $genomeCount,"\t";
		print $hash{$key}{'seq'},"\n";
		}
sub checkClassed{
	my $classMatchFile = "$dir/embl.new_genomes/taxonomyClassification/combined.taxonomy.ClassMatch.txt";
	open (my $fh1, "<","$classMatchFile") or die ("make sure taxonomy matches is in the path $dir/taxonomyClassification/ $!\n");
	my $target = shift; 
	while  (my $line = <$fh1>){
		chomp $line;
		if ($line =~ /^$target/){
			my ($file, $class) = split(" ", $line);
			return "$class";
			}
		}
	}

#print Dumper(\%hash);
#print Dumper(\%seen);
