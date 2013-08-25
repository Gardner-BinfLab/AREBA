#!/usr/bin/python2.7

'''
Created on 12/08/2013

Algorithm:

foreach(region){if(annotated && expression>genomicMedian){$count++}elsif(unannotated && expression<genomicMedian){$count++} $total++;}
concordance=$count/$total;


Multiple GFF files can be used because it is important to mark annotated regions.


@author: suu13
'''
import numpy
import argparse

def GenomeMedian(PlotObject): 
    GenomeMedian=[]
    LenRange=len(PlotObject)
    for i in range(0,LenRange):
        GenomeMedian.append(int(PlotObject[i].split()[0]))
        GenomeMedian.append(int(PlotObject[i].split()[1]))
    return(int(numpy.median(GenomeMedian)))






def GenomeConcordance(PlotObject,GffObject,GenomeM):
    LenRange=len(PlotObject)
    GenomeScaffold=numpy.zeros(shape=(LenRange,2)) #2 sutun ve genome boyu kadar 0 dolu 2 boyutlu array olusturmak
    for GffLine in GffObject: #bu kisim genomdaki her bolgeye annotation varsa 1 yoksa 0 koyuyor
        try: #basinda info kismi olan GFF dosyalari icin ise yariyor 
            if GffLine.split()[6]=='-':
                for i in range(int(GffLine.split()[3])-1,int(GffLine.split()[4])): #plot dosyalari sifirdan basliyor diye -1 yaptim
                    GenomeScaffold[i,0]=1
            else:
                for i in range(int(GffLine.split()[3])-1,int(GffLine.split()[4])):
                    GenomeScaffold[i,1]=1                 
        except:
            pass
    
    concordance=0
    for i in range(0,LenRange): #genome medianindan buyukse ve annotation varsa
        if (GenomeScaffold[i,0]==1 and (GenomeM < int(PlotObject[i].split()[0]))) or (GenomeScaffold[i,0]==0 and (GenomeM > int(PlotObject[i].split()[0]))):
            concordance += 1
        if (GenomeScaffold[i,1]==1 and (GenomeM < int(PlotObject[i].split()[1]))) or (GenomeScaffold[i,1]==0 and (GenomeM > int(PlotObject[i].split()[1]))):
            concordance += 1
    return (float(concordance)/float(2*LenRange))


def main():
    try:
        with open(args.plot) as PlotFile:
            PlotObject=PlotFile.readlines()
        GffObject=[]
        for gff in args.gff:
            with open(gff) as GffFile:
                GffObject=GffObject + GffFile.readlines() #concatenate input gff files 
    except:
        print "File Read Error"
        return
    
    GenomeM=GenomeMedian(PlotObject)
    print "Concordance: " + str(GenomeConcordance(PlotObject,GffObject,GenomeM))









if __name__ == '__main__':
    #parse command line arguments
    Argument_Parser=argparse.ArgumentParser(prog="Concordance.py")
    Argument_Parser.add_argument('-gff',type=str,help="GFF files of annotations",required=True,action="append")
    Argument_Parser.add_argument('-plot',type=str,help="Transcriptome plot file",required=True)
    args=Argument_Parser.parse_args()
    main()