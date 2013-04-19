#!/usr/bin/perl 

#generate_max_plot.pl: For each input *.plot file generate a max.plot
#                      file. For each nucleotide position (on each
#                      strand if data is stranded) record the maximum
#                      observed read-depth.
#                      

use warnings;
use strict;
use Getopt::Long;
use IO::File;
use List::Util qw< min max >;

my ($help, $emblFile, @gffs, @plotFiles, $verbose);

&GetOptions( 
    "h|help"              => \$help,
    "e|embl=s"            => \$emblFile,
    "p|plot=s@"           => \@plotFiles,
    "v|verbose"           => \$verbose
    );

if( $help ) {
    &help();
    exit(1);
}
elsif (not -s $emblFile){
    print "FATAL: no EMBL file given\n";
    &help();
    exit(1);
}

#0. Validate data:
my $emblLength = 0; 
$emblLength = lengthEmbl($emblFile) if defined($emblFile);

my ($plotLengthPrev,$plotPrev,$plotLength);
foreach my $pf (@plotFiles){
    $plotLength = lengthPlot($pf);
    print "WARNING: [$pf] and [$emblFile] are different lengths [$plotLength],[$emblLength]\n" if ($emblLength>0 && $emblLength!=$plotLength);
    print "WARNING: [$pf] and [$plotPrev] are different lengths [$plotLength],[$plotLengthPrev]\n" if (defined($plotLengthPrev) && $plotLengthPrev!=$plotLength);
    
    $plotLengthPrev=$plotLength;
    $plotPrev=$pf;
}

exit(0);

#1. 
my @fileHandles; 
my $i=0;
foreach my $pf (@plotFiles){
    $fileHandles[$i] = IO::File->new();
    $fileHandles[$i] -> open("< $pf") or die "FATAL: failed to plotfile [$pf]\n[$!]";     
    $i++;
}


for ($i=0; $i<$plotLength; $i++){
    my @max=();
    foreach my $fh (@fileHandles){
	
	my $line=<$fh>;
	my @line=split(/\s+/,$line); 
	@max=@line if (not @max);
	for (my $j=0; $j<scalar(@line); $j++){
	    $max[$j]=max($max[$j], $line[$j]); 
	}
	
    }
    my $printL=join("\t", @max);
    print $printL . "\n";
    
}

exit(0);

######################################################################
sub lengthPlot {

    my $plotFile=shift; 
    
    my $fhP = IO::File->new();
    $fhP -> open("wc -l $plotFile |");
    my $length=0;
    while(my $ss = <$fhP>){
	if($ss=~/(\d+)\s+\S+/){
	    $length=$1;
	}
    }
    
    print "WARNING: [wc -l $plotFile] returns an apparent length of [$length]\n" if ($length==0);    
    return $length;    
}

sub lengthEmbl {
    
    my $emblF=shift; 
    
    my $fhE = IO::File->new();
    $fhE -> open("esl-seqstat $emblF |");
    my $lengthE=0;
    while(my $ss = <$fhE>){
	if($ss=~/Total # residues:\s+(\d+)/){
	    $lengthE=$1;
	}
    }
    
    print "WARNING: [esl-seqstat $emblF] returns an apparent length of [$lengthE]\n" if ($lengthE==0);    
    return $lengthE;     
}


sub help {
    print STDERR <<EOF;

generate_max_plot.pl: for an input embl file and a "plot" files
                      generate a maximum read-depth plot file.

Usage:   find_noncoding_peaks.pl -e <embl> -p <plot1> -p <plot2> -p <plot3> ... 
Options:       -h|--help                     Show this help.

               -e|--embl <embl>     EMBL file [optional, used to validate plot files are the correct length]
	       -p|--plot <plot>     Plot file, read-depths per genome position, formatted for Artemis. For multiple plot files use additional -p\'s.
	       -v|--verbose         Print LOTS of stuff [debug option]

EOF
}


