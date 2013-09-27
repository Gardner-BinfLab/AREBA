#!/usr/bin/python2.7

'''
Created on 12/08/2013

Concordance V2 script to draw ROC plots

Revised 20/09/2013


NOTE:Multiple GFF files can be used because it is important to mark annotated regions.


@author: suu13
'''
import numpy
import argparse
from collections import OrderedDict
from joblib import Parallel,delayed

def Expression_Array(PlotObject):
    EA=[]
    for i in range(0,len(PlotObject)):
        EA.append(int(PlotObject[i].split()[0])+int(PlotObject[i].split()[1]))
    EA.sort(reverse=True)
    
    return list(OrderedDict.fromkeys(EA)) #return ordered expression list to test each value
    
def ROC_Table(PlotObject,GffObject):
    LenRange=len(PlotObject)
    GenomeScaffold=Genome_Annotation_Marker(GffObject,LenRange,PlotObject)
    ExpressionArray=Expression_Array(PlotObject)
    print "Expression\tTruePositive\tFalsePositive\tTrueNegative\tFalseNegative"
    
    
    Parallel(n_jobs=12,verbose=5)(delayed(GenomeConcordanceStats)(i,GenomeScaffold,LenRange) for i in ExpressionArray) #12 process run in parallel
    #for i in ExpressionArray:
    #    GenomeCS=GenomeConcordanceStats(i,GenomeScaffold,LenRange)
    #    print "%d\t%d\t%d\t%d\t%d" % (i,GenomeCS[0],GenomeCS[1],GenomeCS[2],GenomeCS[3]) 
        
        
    

def Genome_Annotation_Marker(GffObject,LenRange,PlotObject):
    GenomeScaffold=numpy.zeros(shape=(LenRange,2))
    for GffLine in GffObject:
        try:                 
            for i in range(int(GffLine.split()[3])-1,int(GffLine.split()[4])):
                GenomeScaffold[i,0]=1
        except:
            pass
    for i in range(0,LenRange):
        GenomeScaffold[i,1]=int(PlotObject[i].split()[0])+int(PlotObject[i].split()[1])
     
    return GenomeScaffold #return a genome scaffold with annotations marked, assume unstranded

def GenomeConcordanceStats(ExpressionTreshold,GenomeScaffold,LenRange): #second function to calculate concordance
    TruePositive=0
    FalsePositive=0
    TrueNegative=0
    FalseNegative=0
    for i in range(0,LenRange):
        if (GenomeScaffold[i,0]==1) and (ExpressionTreshold <= GenomeScaffold[i,1]): #assume unstranded plot
            TruePositive +=1
        elif (GenomeScaffold[i,0]==0) and (ExpressionTreshold < GenomeScaffold[i,1]):
            FalsePositive +=1
        elif (GenomeScaffold[i,0]==0) and (ExpressionTreshold > GenomeScaffold[i,1]):
            TrueNegative +=1
        else:
            FalseNegative +=1
        #if (GenomeScaffold[i,0]==1 or GenomeScaffold[i,1]==1) and (ExpressionTreshold >= ExpressionOnLocation):
        #    FalseNegative +=1
    
    #PPV=float(TruePositive/float(TruePositive+FalsePositive+1))
    #Specificity=float(TrueNegative/float(TrueNegative+FalsePositive+1))
    
    
    output="%d\t%d\t%d\t%d\t%d" % (ExpressionTreshold,TruePositive,FalsePositive,TrueNegative,FalseNegative)
    with open(str(ExpressionTreshold)+".exp","w") as tresh:
        tresh.write(output)
    return



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
    
    ROC_Table(PlotObject,GffObject)
    
    
    
if __name__ == '__main__':
    #parse command line arguments
    Argument_Parser=argparse.ArgumentParser(prog="Concordance.py")
    Argument_Parser.add_argument('-gff',type=str,help="GFF files of annotations",required=True,action="append")
    Argument_Parser.add_argument('-plot',type=str,help="Transcriptome plot file",required=True)
    args=Argument_Parser.parse_args()
    main()