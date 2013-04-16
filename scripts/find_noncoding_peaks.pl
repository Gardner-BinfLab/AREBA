#!/usr/bin/perl 

#find_noncoding_peaks.pl: for an input embl file, additional anotations and
#                         a "plot" file, 

use warnings;
use strict;
use Getopt::Long;
use IO::File;
use List::Util qw< min max >;

my ($help, $emblFile, @gffs, $plotFile, $verbose);

&GetOptions( 
    "h|help"              => \$help,
    "e|embl=s"            => \$emblFile,
    "g|gff=s@"            => \@gffs,
    "p|plot=s"            => \$plotFile,
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

#0. Initialize
my $length = lengthEmbl($emblFile);
print "Genome length: [$length]\n" if defined($verbose);

#1. generate GFF from EMBL file
print "Running embl2gff on [$emblFile]\n" if defined($verbose);
my $emblGff = embl2gff($emblFile); 
print "Finished embl2gff, adding [$emblGff] to \@gffs\n" if defined($verbose);
push(@gffs, $emblGff);

#2. mask annotated portions of the genome
#MOVE TO A FUNCTION:
my $gffsS = join(' ',    @gffs); 

my $fhIn  = IO::File->new();
my $fhOut = IO::File->new();

print "Running [sortGffs.pl $gffsS > $$\-sorted.gff && collapseGffs.pl -g $$\-sorted.gff |]\n" if defined($verbose); 

$fhIn -> open( "sortGffs.pl $gffsS > $$\-sorted.gff && collapseGffs.pl -g $$\-sorted.gff |" ) or die "FATAL: failed to run [sortGffs.pl $gffsS > $$\-sorted.gff && collapseGffs.pl -g $$\-sorted.gff |]";
my $annMask="annotation_mask.txt";
open($fhOut, "> $annMask" );
my $posn=1; 
while(my $cur = <$fhIn>) {
    next if ($cur=~/^#/);
    my @cur = split(/\t/, $cur);
    next if (scalar(@cur)!=9);
    for (;$posn<$cur[3]; $posn++){
	print $fhOut "0\n";
    }
    
    for (;$posn<$cur[4]-1; $posn++){
	print $fhOut "1\n";
    }
}

for (;$posn<$length+1; $posn++){
    print $fhOut "0\n";
}

$fhIn ->close();
$fhOut->close();

unlink($emblGff);
unlink("$$\-sorted.gff");

#3. print a ranked list of unannotated regions based upon expression level (max, mean, median)
my @expressionVals=();
my $unannExp = "$plotFile\-unannotated-expression.txt";
open(PT, "< $plotFile") or die "FATAL: failed to open [$plotFile]!\n[$!]" ;
open(AM, "< $annMask")  or die "FATAL: failed to open [$annMask]!\n[$!]" ;
open(UT, "> $unannExp") or die "FATAL: failed to open [$unannExp]!\n[$!]" ;

printf UT "#start\tend\tlength\tmean\tmedian\tmax\n"; 	
my ($start,$sum, $i)=(0,0,0);

while(defined(my $mask=<AM>) && defined(my $plotVals=<PT>)){   
    chomp($mask);
    chomp($plotVals);
    if($mask==0){
	$start=$i if ($start==0);
	my @plotValsI=split(/\s+/, $plotVals);	
	push(@expressionVals, max(@plotValsI));
	$sum+=max(@plotValsI);
    }
    elsif($mask==1 && scalar(@expressionVals)){
	my $mean   = $sum/($i-$start+1);
	my $median = median(\@expressionVals);
	my $max    = max(@expressionVals);
	printf UT "$start\t$i\t%d\t%0.2f\t$median\t$max\n", $i-$start+1, $mean;
	@expressionVals=();
	($start,$sum)=(0,0);
    }
    $i++;
}

if(scalar(@expressionVals)){    
    my $mean  = $sum/($i-$start+1);
    my $median= median(\@expressionVals);
    my $max = max @expressionVals;	
    printf UT "$start\t$i\t%d\t%0.2f\t$median\t$max\n", $i-$start+1, $mean;
    @expressionVals=();
    ($start,$sum)=(0,0);
}

close(AM);
close(PT);
close(UT);

exit(0);

######################################################################
#median
sub median {
    my $vals=shift;
    my @sorted_values = sort by_number @{$vals};
    my $mid = int @sorted_values/2;
    my $median;
    if (@sorted_values % 2) {
	$median = $sorted_values[ $mid ];
    } else {
	$median = ($sorted_values[$mid-1] + $sorted_values[$mid])/2;
    } 
    
    return $median;
}

######################################################################
#by_number
sub by_number {
    if ($a < $b) {
        return -1;
    } 

elsif ($a == $b) {
        return 0;
    } elsif ($a > $b) {
        return 1;
    }
}

######################################################################
#embl2gff: extract annotations from an EMBL file, print them in GFF
sub embl2gff {
    
    my $emblFile=shift;
    
    print "Opening [$emblFile]\n" if defined($verbose);
    open(EM, "< $emblFile") or die "FATAL: failed to open [$emblFile]!\n[$!]" ;
    my(%store);#$type,$s,$e,$strand,$gene,$product,$notes);
    my $seqId = 'undefined';
    
    
    
    my $gffOutFile="$$.gff";
    my $fh = IO::File->new();
    open($fh, ">", $gffOutFile );
    
    print "Printing EMBL GFF to [$gffOutFile]\n" if defined($verbose); 
    
    while(my $l=<EM>){
	
	if($l=~/^ID\s+(\S+;)/){
	    $seqId=$1;
	}
	elsif($l=~/^FT\s{3}(\S+)\s+(\d+)\.\.(\d+)/){
	    next if $1 eq 'source';
	    printGff($seqId,\%store,$fh) if (defined $seqId && defined $store{'s'} && defined $store{'e'} && defined $store{'strand'}); 
	    undef %store;
	    ($store{'type'},$store{'s'},$store{'e'},$store{'strand'})=($1,$2,$3,'+');
	}
	elsif($l=~/^FT\s{3}(\S+)\s+complement\((\d+)\.\.(\d+)\)/){
	    printGff($seqId,\%store,$fh) if (defined $seqId && defined $store{'s'} && defined $store{'e'} && defined $store{'strand'}); 
	    undef %store;
	    ($store{'type'},$store{'s'},$store{'e'},$store{'strand'})=($1,$2,$3,'-');
	}
	elsif($l=~/^FT\s+\/gene="(.*)/){
	    $store{'gene'}=$1;
	    $store{'gene'}=~s/\"$//;
	}
	elsif($l=~/^FT\s+\/product="(.*)/){
	    $store{'product'}=$1;
	    $store{'product'}=~s/\"$//;
	}
	elsif($l=~/^FT\s+\/note="(.*)/){
	$store{'note'}=$1;
	$store{'note'}=~s/\"$//;
	}
	elsif($l=~/^FT\s+\/locus_tag="(.*)/){
	    $store{'locus_tag'}=$1;
	    $store{'locus_tag'}=~s/\"$//;
	}
	
    }
    close(EM);
    
    $fh->close;
    
    return $gffOutFile;
    
    
}

sub printGff {
    my ($seqId,$store,$fh)=@_;
    $store->{'type'}='gene' if (not defined $store->{'type'});
    my $phase='.';
    $phase=0 if $store->{'type'}=~/CDS/;
    my ($attributes,$spacer)=('','');
    $attributes .= 'Name='    . $store->{'gene'}                     if (defined $store->{'gene'});
    $spacer = ';' if length($attributes)>0;
    $attributes .= $spacer . 'Product="' . $store->{'product'} . '"' if (defined $store->{'product'});
    $spacer = ';' if length($attributes)>0;
    $attributes .= $spacer . 'Note="'    . $store->{'notes'}   . '"' if (defined $store->{'notes'});
    $attributes .= $spacer . 'Locus_Tag="'    . $store->{'locus_tag'}   . '"' if (defined $store->{'locus_tag'});
    $attributes .= ';';

    $fh -> print("$seqId\tEMBL\t".$store->{'type'}."\t".$store->{'s'}."\t".$store->{'e'}."\t.\t".$store->{'strand'}."\t.\t$attributes\n");
    #print "gff:[$seqId\tEMBL\t".$store->{'type'}."\t".$store->{'s'}."\t".$store->{'e'}."\t.\t".$store->{'strand'}."\t.\t$attributes]\n" if defined($verbose);
    return 0;
}

sub lengthEmbl {
    
    my $emblFile=shift; 
    
    my $fhE = IO::File->new();
    $fhE -> open("esl-seqstat $emblFile |");
    my $length=0;
    while(my $ss = <$fhE>){
	if($ss=~/Total # residues:\s+(\d+)/){
	    $length=$1;
	}
    }
    
    print "WARNING: [esl-seqstat $emblFile] returns an apparent length of [$length]\n" if ($length==0);    
    return $length;     
}

sub help {
    print STDERR <<EOF;

find_noncoding_peaks.pl: for an input embl file, additional anotations and
                         a "plot" file,

Usage:   find_noncoding_peaks.pl -e <embl> -p <plot> -g <gff1> -g <gff2> -g <gff3> ... 
Options:       -h|--help                     Show this help.

               -e|--embl <embl>     EMBL file
	       -g|--gff  <gff>      Give GFF file names as input. For multiple gff files use additional -g\'s.
	       -p|--plot <plot>     Plot file, read-depths per genome position, formatted for Artemis. 
	       -v|--verbose         Print LOTS of stuff [debug option]

EOF
}


