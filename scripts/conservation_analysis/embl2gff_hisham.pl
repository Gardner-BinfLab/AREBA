use strict;
use warnings;
use File::Find;
use Bio::SeqIO;
use Data::Dumper;


#A program to convert from embl format into GFF format using bioPerl.
#Can be useful for many other projects too.
#Syntax:  perl embl2gff_hisham.pl | grep RUF > embl.new_genomes/rufEntries.gff 
#biohisham@gmail.com

my @files; #embl files
find(\&wanted, "./embl.new_genomes");
sub wanted{
	if ($File::Find::name =~ /\.\/(.*\.embl)/){
		push @files, $1;
		};	
	}

foreach my $element (@files){ #initialize Bio::SeqIO objects and convert each file to its equivalent gff counterpart.
	#print $element, $/;
	open (my $gff_FH,">","$element.gff") or die("could not open file $!\n");
	my $in = Bio::SeqIO->new(-file => $element, -format => 'embl');
	while(my $seq = $in->next_seq){
		for my $feat ($seq->top_SeqFeatures){
			print $gff_FH $feat->gff_string, "\n"; #activate when printing to files
			print STDOUT $feat->gff_string, "\n"; # #activate when grepping and dumping to rufEntries.gff
			}
		}
	}

#Clean the entries and retain the first, fourth and fifth column from the gff file
# perl -lane 'print "$F[0]\t$F[3]\t$F[4]"' embl.new_genomes/rufEntries.gff >embl.new_genomes/rufEntries_clean.gff

#run extractSeqs.pl
