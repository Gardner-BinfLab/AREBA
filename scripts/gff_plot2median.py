#!/usr/bin/env python2.7
'''
Created on 9/05/2013

@author: suu13
'''

import numpy
import argparse
import random

def GenomeConcordance(PlotObject,GffObject,GenomeM):
    LenRange=len(PlotObject)
    GenomeScaffold=numpy.zeros(shape=(LenRange,2)) #2 sutun ve genome boyu kadar 0 dolu 2 boyutlu array olusturmak
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
    
    concordance=0
    for i in range(0,LenRange):
        if (GenomeScaffold[i,0]==1 and GenomeM < int(PlotObject[i].split()[0])) or (GenomeScaffold[i,0]==0 and GenomeM > int(PlotObject[i].split()[0])):
            concordance += 1
        if (GenomeScaffold[i,1]==1 and GenomeM < int(PlotObject[i].split()[1])) or (GenomeScaffold[i,1]==0 and GenomeM > int(PlotObject[i].split()[1])):
            concordance += 1
    
    return (concordance/2*LenRange)
    
    

def Rfam_Pfam_Stats(PlotObject,GffObject): #to calculate reads per Rfam and Pfam annotation
    LenRange=len(PlotObject)
    GenomeScaffold=numpy.zeros(shape=(LenRange,1),dtype=int) #1 sutun ve genome boyu kadar 0 dolu 1 boyutlu array olusturmak
    for GffLine in GffObject:
        try: #basinda info kismi olan GFF dosyalari icin ise yariyor 
            if GffLine.split()[1]=='Rfam':
                for i in range(int(GffLine.split()[3])-1,int(GffLine.split()[4])): #plot dosyalari sifirdan basliyor diye -1 yaptim
                    GenomeScaffold[i,0]=1 #Rfam icin 1
            elif GffLine.split()[1]=='Pfam':
                for i in range(int(GffLine.split()[3])-1,int(GffLine.split()[4])):
                    GenomeScaffold[i,0]=2 #Pfam icin 2
            else:
                for i in range(int(GffLine.split()[3])-1,int(GffLine.split()[4])):
                    GenomeScaffold[i,0]=0 # geri kalanlar 0              
        except:
            pass
    Rfam_Reads=0    
    Pfam_Reads=0
    Total_Reads=0
    for i in range(0,LenRange):
        Total_Reads=Total_Reads+int(PlotObject[i].split()[0])+int(PlotObject[i].split()[1])
        if (GenomeScaffold[i,0]==1):
            Rfam_Reads=Rfam_Reads+int(PlotObject[i].split()[0])+int(PlotObject[i].split()[1])
        elif (GenomeScaffold[i,0]==2):
            Pfam_Reads=Pfam_Reads+int(PlotObject[i].split()[0])+int(PlotObject[i].split()[1])
        
    return (Rfam_Reads,Pfam_Reads,Total_Reads)
        




def GenomeMedian(PlotObject): 
    GenomeMedian=[]
    LenRange=len(PlotObject)
    for i in range(0,LenRange):
        GenomeMedian.append(int(PlotObject[i].split()[0])+int(PlotObject[i].split()[1]))
    return(int(numpy.median(GenomeMedian))) #find the genomic median, assume unstranded
        
    

#function to calculate normalized Median of strands, assumes both stranded or non-stranded plot files, uses median tuple
def MedianNormalize(GffLine,PlotObject,MedianTup):
    LenRange=len(PlotObject)
    TranscriptLen=abs(int(GffLine.split()[3])-int(GffLine.split()[4]))
    MedianMedArrayStrand=[]
    MedianMedArrayCombined=[]
    for _ in range(1,201): #number of randomization
        try:
            MedArrayStrand=[]
            MedArrayCombined=[]
            RandomStartLoc=random.randint(1,LenRange) #random start location selection
            RandomStrand=random.randint(0,1) #random strand selection
            if RandomStrand==0:
                for i in range(RandomStartLoc,RandomStartLoc+TranscriptLen): #plot dosyalari sifirdan basliyor diye -1 yaptim
                    MedArrayStrand.append(int(PlotObject[i].split()[0]))
                    MedArrayCombined.append(int(PlotObject[i].split()[0])+int(PlotObject[i].split()[1]))
            else:
                for i in range(RandomStartLoc,RandomStartLoc+TranscriptLen):
                    MedArrayStrand.append(int(PlotObject[i].split()[1]))
                    MedArrayCombined.append(int(PlotObject[i].split()[0])+int(PlotObject[i].split()[1]))
            MedianMedArrayStrand.append(int(numpy.median(MedArrayStrand)))
            MedianMedArrayCombined.append(int(numpy.median(MedArrayCombined)))
        except:
            pass
    return (MedianTup[0]/(1+int(numpy.median(MedianMedArrayStrand))),MedianTup[1]/(1+int(numpy.median(MedianMedArrayCombined))))
            



def MedianTranscript(GffLine,PlotObject): #function to calculate Median of strands, assumes both stranded or non-stranded plot files
    MedArrayStrand=[]
    MedArrayCombined=[]
    if GffLine.split()[6]=='-':
        for i in range(int(GffLine.split()[3])-1,int(GffLine.split()[4])): #plot dosyalari sifirdan basliyor diye -1 yaptim
            MedArrayStrand.append(int(PlotObject[i].split()[1])) #AREBA da hata varmis, duzelttim plot dosyasinda forward ilk, reverse ikinci sutun
            MedArrayCombined.append(int(PlotObject[i].split()[0])+int(PlotObject[i].split()[1]))
    else:
        for i in range(int(GffLine.split()[3])-1,int(GffLine.split()[4])):
            MedArrayStrand.append(int(PlotObject[i].split()[0])) #AREBA da hata varmis, duzelttim plot dosyasinda forward ilk, reverse ikinci sutun
            MedArrayCombined.append(int(PlotObject[i].split()[0])+int(PlotObject[i].split()[1]))
    
    return (int(numpy.median(MedArrayStrand)),int(numpy.median(MedArrayCombined)),int(numpy.mean(MedArrayStrand)),int(numpy.mean(MedArrayCombined)),int(numpy.max(MedArrayStrand)),
            int(numpy.max(MedArrayCombined)),int(numpy.min(MedArrayStrand)),int(numpy.min(MedArrayCombined)))

def Median_Window_Transcript(Window,PlotObject): #function to calculate Median of strands, assumes both stranded or non-stranded plot files
    MedArrayCombined=[]
    if(Window[0]<Window[1]):
        w=range(Window[0]-1,Window[1])
    else:
        w=range(Window[1]-1,Window[0])
    
    
    for i in w: #plot dosyalari sifirdan basliyor diye -1 yaptim
        MedArrayCombined.append(int(PlotObject[i].split()[0])+int(PlotObject[i].split()[1]))
    return (int(numpy.median(MedArrayCombined)))





def main():
#   try:
    with open(args.plot) as PlotFile:
        PlotObject=PlotFile.readlines()
    
    if(args.mediangenome==True):
        print "Median of %s:%d" %(args.plot,GenomeMedian(PlotObject))
        return
    elif(args.window!=None):
        print Median_Window_Transcript(args.window,PlotObject)
        return
    
    GffObject=[]
    for gff in args.gff:
        with open(gff) as GffFile:
            GffObject=GffObject + GffFile.readlines() #concatenate input gff files 
#    except:
#       print "File Read Error"
#        return
    
   
    
    if(args.rfampfamstats==True):
        Rfam_Pfam_S=(Rfam_Pfam_Stats(PlotObject,GffObject))
        print "%i\t%i\t%i" % (Rfam_Pfam_S)
    else:
        for GffLine in GffObject:
            #try: #basinda info kismi olan GFF dosyalari icin ise yariyor
            MedianTup=MedianTranscript(GffLine,PlotObject)
            #print "%s;Median_depth_one_strand:$%d$,two_strands_combined:$%d" %(GffLine.strip(),MedianTup[0],MedianTup[1])
            print "%s;Median_depth_one_strand:%d,two_strands_combined:%d,Mean_depth_one_strand:%d,two_strands_combined:%d,Max_strand:%d,combined:%d,Min_strand:%d,combined:%d" %(
                    GffLine.strip(),MedianTup[0],MedianTup[1],MedianTup[2],MedianTup[3],MedianTup[4],MedianTup[5],MedianTup[6],MedianTup[7])
            #except:
            #    pass

    
    """    
    else:
        GenomeM=GenomeMedian(PlotObject)
        for GffLine in GffObject:
            try: #basinda info kismi olan GFF dosyalari icin ise yariyor
                MedianTup=MedianTranscript(GffLine,PlotObject)
                if (args.additional==False):
                    print "%s;Median_depth_stranded:$%d$,combined:$%d$,(Genome_Median:$%d$)" %(GffLine.strip(),MedianTup[0],MedianTup[1],GenomeM)
                else:
                    MedianNorm=MedianNormalize(GffLine,PlotObject,MedianTup)
                    print "%s;Median_depth_stranded:$%d$,combined:$%d$,Normalized:$%d$%d$" %(GffLine.strip(),MedianTup[0],MedianTup[1],MedianNorm[0],MedianNorm[1])
    
            except:
                pass
    """ 

            
    return
            
if __name__ == '__main__':
    #parse command line arguments
    Argument_Parser=argparse.ArgumentParser(prog="gff_plot2median.py")
    Argument_Parser.add_argument('-gff',type=str,help="GFF files of annotations",action="append")
    Argument_Parser.add_argument('-plot',type=str,help="Transcriptome plot file",required=True)
    Argument_Parser.add_argument('-additional',action='store_true',help="Additional")#to add additional statistics
    Argument_Parser.add_argument('-mediangenome',action='store_true',help="Return median of genome only")
    Argument_Parser.add_argument('-rfampfamstats',action='store_true',help="Rfam Pfam Stats")
    Argument_Parser.add_argument('-window',type=int,nargs=2,help="Window to extract from sequences")
    args=Argument_Parser.parse_args()
    main()

        