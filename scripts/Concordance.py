#!/usr/bin/python2.7

'''
Created on 12/08/2013

Old Algorithm:

foreach(region){if(annotated && expression>genomicMedian){$count++}elsif(unannotated && expression<genomicMedian){$count++} $total++;}
concordance=$count/$total;



Revised 20/09/2013


NOTE:Multiple GFF files can be used because it is important to mark annotated regions.


@author: suu13
'''
import numpy
import argparse

 
        
    

def GenomeMedian(PlotObject): 
    GenomeMedian=[]
    LenRange=len(PlotObject)
    for i in range(0,LenRange):
        GenomeMedian.append(int(PlotObject[i].split()[0])+int(PlotObject[i].split()[1]))
    return(int(numpy.median(GenomeMedian))) #find the genomic median, assume unstranded

def Genome_Annotation_Marker(GffObject,PlotObject):
    GenomeScaffold=numpy.zeros(shape=(len(PlotObject),2))
    for GffLine in GffObject:
        try:                 
            for i in range(int(GffLine.split()[3])-1,int(GffLine.split()[4])):
                GenomeScaffold[i,0]=1
        except:
            pass
    for i in range(0,len(PlotObject)):
        GenomeScaffold[i,1]=int(PlotObject[i].split()[0])+int(PlotObject[i].split()[1])
     
    return GenomeScaffold #return a genome scaffold with annotations marked, and expression summed: assume unstranded data

"""
def Genome_Annotation_Marker(GffObject,PlotLength):
    GenomeScaffold=numpy.zeros(shape=(PlotLength,2))
    for GffLine in GffObject:
        try: #basinda info kismi olan GFF dosyalari icin ise yariyor 
            if GffLine.split()[6]=='-':
                for i in range(int(GffLine.split()[3])-1,int(GffLine.split()[4])): #plot dosyalari sifirdan basliyor diye -1 yaptim
                    GenomeScaffold[i,0]=1
            else:
                for i in range(int(GffLine.split()[3])-1,int(GffLine.split()[4])):
                    GenomeScaffold[i,1]=1                   
        except:
            pass 
    return GenomeScaffold
 

def GenomeConcordanceStats(PlotObject,GffObject,ExpressionTreshold,GenomeScaffold): #second function to calculate concordance
    LenRange=len(PlotObject)
    TruePositive=0
    FalsePositive=0
    TrueNegative=0
    FalseNegative=0
    for i in range(0,LenRange):
        ExpressionOnLocation=int(PlotObject[i].split()[0])+int(PlotObject[i].split()[1])
        if (GenomeScaffold[i,0]==1 or GenomeScaffold[i,1]==1) and (ExpressionTreshold < ExpressionOnLocation): #assume unstranded plot
            TruePositive +=1
        elif (GenomeScaffold[i,0]==0 and GenomeScaffold[i,1]==0) and (ExpressionTreshold < ExpressionOnLocation):
            FalsePositive +=1
        elif (GenomeScaffold[i,0]==0 and GenomeScaffold[i,1]==0) and (ExpressionTreshold >= ExpressionOnLocation):
            TrueNegative +=1
        else:
            FalseNegative +=1
        #if (GenomeScaffold[i,0]==1 or GenomeScaffold[i,1]==1) and (ExpressionTreshold >= ExpressionOnLocation):
        #    FalseNegative +=1
    
    PPV=float(TruePositive/float(TruePositive+FalsePositive+1))
    Specificity=float(TrueNegative/float(TrueNegative+FalsePositive+1))
    
    return [PPV,Specificity]
"""

def GenomeConcordanceStats(ExpressionTreshold,GenomeScaffold): #second function to calculate concordance
    TruePositive=0
    FalsePositive=0
    TrueNegative=0
    FalseNegative=0
    for i in range(0,len(GenomeScaffold)):
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
    
    PPV=float(TruePositive/float(TruePositive+FalsePositive))
    Specificity=float(TrueNegative/float(TrueNegative+FalsePositive))
    
    return [PPV,Specificity]

def GenomeConcordance(PlotObject,GffObject,GenomeM): #old function to calculate concordance
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
    
    if(args.roc==False):
        GenomeM=GenomeMedian(PlotObject)
        GenomeScaffold=Genome_Annotation_Marker(GffObject,PlotObject)
        print "Concordance: " + str(GenomeConcordance(PlotObject,GffObject,GenomeM))

    else:
        GenomeM=GenomeMedian(PlotObject)
        GenomeScaffold=Genome_Annotation_Marker(GffObject,PlotObject)                
        GenomeCS=GenomeConcordanceStats(GenomeM,GenomeScaffold)
        print "PPV & Specificity:%f\t%f" % (GenomeCS[0],GenomeCS[1])
    
    
    









if __name__ == '__main__':
    #parse command line arguments
    Argument_Parser=argparse.ArgumentParser(prog="Concordance.py")
    Argument_Parser.add_argument('-gff',type=str,help="GFF files of annotations",required=True,action="append")
    Argument_Parser.add_argument('-plot',type=str,help="Transcriptome plot file",required=True)
    Argument_Parser.add_argument('-roc',action='store_true',help="ROC_Table")
    args=Argument_Parser.parse_args()
    main()