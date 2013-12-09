#!/bin/bash

#$ -S /bin/bash
#$ -cwd
#cd ${PBS_O_WORKDIR}; 
#if [ `grep "GS" $1 -c` -ge 5000 ]; then esl-weight -f --idf 0.3 $1 > $1.weight; fi;
#if [ `grep "GS" $1 -c` -le 5000 ] &&  [ `grep "GS" $1 -c` -ge 200 ]; then esl-weight -f --idf 0.55 $1 > $1.weight; fi;
if [ ! -f $1.clustal ]; then /home/suu13/bin/biof_converter.py -f $1 -i stockholm -o clustal > $1.clustal; fi;

#if [ ! -f $1*RNAz ]; then /export/home/ppg15/inst/RNAz-2.1/rnaz/RNAz -b $1.weight.clustal 1> $1.weight.clustal.RNAz 2> $1.weight.clustal.RNAz.err; fi;
if [ ! -f $1*RNAz ]; then /home/suu13/progs/RNAz/bin/RNAz -b $1.clustal 1> $1.clustal.RNAz 2> $1.clustal.RNAz.err; fi;

#if [ ! -f $1*RNAcode ]; then RNAcode --tabular --best-only $1.weight.clustal 1> $1.weight.clustal.RNAcode 2> $1.weight.clustal.RNAcode.err; fi;
#if [ ! -f $1*RNAcode ]; then RNAcode --tabular --best-only $1.clustal 1> $1clustal.RNAcode 2> $1.clustal.RNAcode.err; fi;









#biof_converter.py -f $1 -i stockholm -o fasta > $1.fasta;
#cat $1.fasta | awk '{if(/>/) print $1; else print $0;}' | sponge $1.fasta;




