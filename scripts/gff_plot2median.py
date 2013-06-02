#!/usr/bin/python2.7
'''
Created on 9/05/2013

@author: suu13
'''

import numpy
import argparse
import random
from scipy.stats import pearsonr

"""
def PlotStat(AdditionalStat):
    Cumulative=[] #combine stranded median counts of annotations
    for i in range(0,len(AdditionalStat)):
        Cumulative.append(abs(AdditionalStat[i][0]-AdditionalStat[i][1]))
    Cumulative.sort(reverse=True)
    print Cumulative
    return float(sum(Cumulative)/float(len(Cumulative)))
"""        
    

def CorrelationCoef(GffLine,PlotObject): #pearson correlation coef per annotation
    ReverseStrand=[]
    ForwardStrand=[]
    for i in range(int(GffLine.split()[3])-1,int(GffLine.split()[4])):
        ForwardStrand.append(int(PlotObject[i].split()[1]))
        ReverseStrand.append(int(PlotObject[i].split()[0]))
    return pearsonr(ReverseStrand,ForwardStrand)[0]
        
    

def GenomeDepthPerMillion(PlotObject): 
    TotalDepth=0
    #GenomeMedian=[]
    LenRange=len(PlotObject)
    for i in range(0,LenRange):
        TotalDepth=TotalDepth+int(PlotObject[i].split()[0])+int(PlotObject[i].split()[1])
        #GenomeMedian.append(int(PlotObject[i].split()[0]))
        #GenomeMedian.append(int(PlotObject[i].split()[1]))
    return(float(TotalDepth/1000000))
    
    

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
            

#function to calculate Median of strands, assumes both stranded or non-stranded plot files
def MedianTranscript(GffLine,PlotObject):
    MedArrayStrand=[]
    MedArrayCombined=[]
    if GffLine.split()[6]=='-':
        for i in range(int(GffLine.split()[3])-1,int(GffLine.split()[4])): #plot dosyalari sifirdan basliyor diye -1 yaptim
            MedArrayStrand.append(int(PlotObject[i].split()[0]))
            MedArrayCombined.append(int(PlotObject[i].split()[0])+int(PlotObject[i].split()[1]))
    else:
        for i in range(int(GffLine.split()[3])-1,int(GffLine.split()[4])):
            MedArrayStrand.append(int(PlotObject[i].split()[1]))
            MedArrayCombined.append(int(PlotObject[i].split()[0])+int(PlotObject[i].split()[1]))
    
    return (int(numpy.median(MedArrayStrand)),int(numpy.median(MedArrayCombined)))

def main():
    try:
        with open(args.plot) as PlotFile:
            PlotObject=PlotFile.readlines()
        with open(args.gff) as GffFile:
            GffObject=GffFile.readlines()
    except:
        print "File Read Error"
        return
    if(args.additional==True):
        GDPM=GenomeDepthPerMillion(PlotObject)
        AdditionalStat=[] #Additional statistics list to record MedianTuple
    for GffLine in GffObject:
        try: #basinda info kismi olan GFF dosyalari icin ise yariyor
            MedianTup=MedianTranscript(GffLine,PlotObject)
            #NormalizedMedianTup=MedianNormalize(GffLine,PlotObject,MedianTup)
            #print GffLine.strip()+";Median_depth_stranded:"+str(MedianTup[0])+"("+str(float(MedianTup[0]/GDPM))+")"+",combined:"+str(MedianTup[1])+"("+str(float(MedianTup[1]/GDPM))+")"
            #print "%s\tMedian_depth_stranded:%d(DPM:%.2f)(Norm:%d),combined:%d(DPM:%.2f)(Norm:%d),correlation:%.2f" %(GffLine.strip(),MedianTup[0],MedianTup[0]/GDPM,NormalizedMedianTup[0],MedianTup[1],MedianTup[1]/GDPM,NormalizedMedianTup[1],CorCoef)
            if (args.additional==False):
                print "%s;Median_depth_stranded:$%d$,combined:$%d$" %(GffLine.strip(),MedianTup[0],MedianTup[1])
                pass
            else:
                AdditionalStat.append(MedianTup)
                CorCoef=CorrelationCoef(GffLine,PlotObject)
                print "%s;Median_depth_stranded:$%d$(DPM:$%.2f$),combined:$%d$(DPM:$%.2f$),correlation:$%.2f$" %(GffLine.strip(),MedianTup[0],MedianTup[0]/GDPM,MedianTup[1],MedianTup[1]/GDPM,CorCoef)
        except:
            pass
    #if(args.additional==True):
    #    print PlotStat(AdditionalStat)
            
    return
            
if __name__ == '__main__':
    #parse command line arguments
    Argument_Parser=argparse.ArgumentParser(prog="gff_plot2median.py")
    Argument_Parser.add_argument('-gff',type=str,help="GFF file of annotations",required=True)
    Argument_Parser.add_argument('-plot',type=str,help="Transcriptome plot file",required=True)
    Argument_Parser.add_argument('-additional',action='store_true',help="Additional ")#to add additional statistics
    args=Argument_Parser.parse_args()
    main()

        