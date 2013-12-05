#!/usr/bin/python2.7

'''
Created on 12/08/2013

Concordance V2 script to draw ROC plots

Revised 20/09/2013


NOTE:Multiple GFF files can be used because it is important to mark annotated regions.


@author: suu13
'''
from numpy import zeros
from numpy import median
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
        
        
    
def Genome_Annotation_Marker(GffObject,PlotObject):
    GenomeScaffold=zeros(shape=(len(PlotObject),2))
    for GffLine in GffObject:
        try:                 
            for i in range(int(GffLine.split()[3])-1,int(GffLine.split()[4])):
                GenomeScaffold[i,0]=1
        except:
            pass
    for i in range(0,len(PlotObject)):
        GenomeScaffold[i,1]=int(PlotObject[i].split()[0])+int(PlotObject[i].split()[1])
     
    return GenomeScaffold #return a genome scaffold with annotations marked, and expression summed: assume unstranded data

def Genome_Annotation_Marker_Shuffled(GffObject,PlotObject,ShiftSize):
    GenomeScaffold=zeros(shape=(len(PlotObject),2))
    for GffLine in GffObject:
        try:                 
            for i in range(int(GffLine.split()[3])-1,int(GffLine.split()[4])):
                GenomeScaffold[i,0]=1
        except:
            pass
    for i in range(0,len(PlotObject)):
        try:      
            GenomeScaffold[i,1]=int(PlotObject[i+ShiftSize].split()[0])+int(PlotObject[i+ShiftSize].split()[1])
        except:
            GenomeScaffold[i,1]=int(PlotObject[i+ShiftSize-(len(PlotObject)-1)].split()[0])+int(PlotObject[i+ShiftSize-(len(PlotObject)-1)].split()[1])
        
     
    return GenomeScaffold #return a genome scaffold with annotations marked, and expression summed: assume unstranded data

def GenomeConcordanceStats(ExpressionTreshold,GenomeScaffold,LenRange): #second function to calculate concordance
    TruePositive=0
    FalsePositive=0
    TrueNegative=0
    FalseNegative=0
    for i in range(0,LenRange):
        if (GenomeScaffold[i,0]==1) and (ExpressionTreshold <= GenomeScaffold[i,1]): #assume unstranded plot
            TruePositive +=1
        elif (GenomeScaffold[i,0]==0) and (ExpressionTreshold <= GenomeScaffold[i,1]):
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


def GenomeConcordance(GenomeM,GenomeScaffold): #old function to calculate concordance
    concordance=0
    for i in range(0,len(GenomeScaffold)): #genome medianindan buyukse ve annotation varsa
        if (GenomeScaffold[i,0]==1) and (GenomeM <= GenomeScaffold[i,1]):
            concordance += 1
        elif (GenomeScaffold[i,0]==0) and (GenomeM > GenomeScaffold[i,1]):
            concordance += 1
    return (float(concordance)/len(GenomeScaffold))


def GenomeConcordanceShuffled(GenomeM,GenomeScaffold,ShiftSize):
    concordance=0
    for i in range(0,len(GenomeScaffold)): #genome medianindan buyukse ve annotation varsa
        try:
            if (GenomeScaffold[i,0]==1) and (GenomeM <= GenomeScaffold[i+ShiftSize,1]):
                concordance += 1
            elif (GenomeScaffold[i,0]==0) and (GenomeM > GenomeScaffold[i+ShiftSize,1]):
                concordance += 1
        except:
            if (GenomeScaffold[i,0]==1) and (GenomeM <= GenomeScaffold[i+ShiftSize-len(GenomeScaffold),1]):
                concordance += 1
            elif (GenomeScaffold[i,0]==0) and (GenomeM > GenomeScaffold[i+ShiftSize-len(GenomeScaffold),1]):
                concordance += 1
    #print concordance            
    return (float(concordance)/len(GenomeScaffold))

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
    
    #ROC_Table(PlotObject,GffObject)
    
    #print (args.plot)+"\t",
    Concor=[]
    GenomeScaffold=Genome_Annotation_Marker(GffObject,PlotObject)
    Concor.append(GenomeConcordanceShuffled(10,GenomeScaffold,0))
    
    for ShiftSize in range(50000,2000000,50000):
        #GenomeScaffold=Genome_Annotation_Marker_Shuffled(GffObject,PlotObject,ShiftSize)
        Concor.append(str(GenomeConcordanceShuffled(10,GenomeScaffold,ShiftSize)))
    
    print(','.join(map(str,Concor)))
    
if __name__ == '__main__':
    #parse command line arguments
    Argument_Parser=argparse.ArgumentParser(prog="Concordance.py")
    Argument_Parser.add_argument('-gff',type=str,help="GFF files of annotations",required=True,action="append")
    Argument_Parser.add_argument('-plot',type=str,help="Transcriptome plot file",required=True)

    args=Argument_Parser.parse_args()
    main()