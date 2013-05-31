#!/usr/bin/python2.7
'''
Created on 9/05/2013

@author: suu13
'''

import numpy
import argparse
import random

def MedianNormalize(GffLine,PlotObject,MedianTup):
    LenRange=len(PlotObject)
    TranscriptLen=abs(int(GffLine.split()[3])-int(GffLine.split()[4]))
    MedianMeanArrayStrand=[]
    MedianMeanArrayCombined=[]
    for _ in range(1,501):
        try:
            MedArrayStrand=[]
            MedArrayCombined=[]
            RandomStartLoc=random.randint(1,LenRange)
            RandomStrand=random.randint(0,1)
            if RandomStrand==0:
                for i in range(RandomStartLoc,RandomStartLoc+TranscriptLen): #plot dosyalari sifirdan basliyor diye -1 yaptim
                    MedArrayStrand.append(int(PlotObject[i].split()[0]))
                    MedArrayCombined.append(int(PlotObject[i].split()[0])+int(PlotObject[i].split()[1]))
            else:
                for i in range(RandomStartLoc,RandomStartLoc+TranscriptLen):
                    MedArrayStrand.append(int(PlotObject[i].split()[1]))
                    MedArrayCombined.append(int(PlotObject[i].split()[0])+int(PlotObject[i].split()[1]))
            MedianMeanArrayStrand.append(int(numpy.median(MedArrayStrand)))
            MedianMeanArrayCombined.append(int(numpy.median(MedArrayCombined)))
        except:
            pass
    return (MedianTup[0]/(1+int(numpy.median(MedianMeanArrayStrand))),MedianTup[0]/(1+int(numpy.median(MedianMeanArrayCombined))))
            


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
    
    for GffLine in GffObject:
        try: #basinda info kismi olan GFF dosyalari icin ise yariyor
            MedianTup=MedianTranscript(GffLine,PlotObject)
            NormalizedMedianTup=MedianNormalize(GffLine,PlotObject,MedianTup)
            print GffLine.strip()+";Median_depth_stranded:"+str(MedianTup[0])+"("+str(NormalizedMedianTup[0])+")"+",combined:"+str(MedianTup[1])+"("+str(NormalizedMedianTup[1])+")"
        except:
            pass
            
            
if __name__ == '__main__':
    Argument_Parser=argparse.ArgumentParser(prog="gff_plot2median.py")
    Argument_Parser.add_argument('-gff',type=str,help="GFF file of annotations",required=True)
    Argument_Parser.add_argument('-plot',type=str,help="Transcriptome plot file",required=True)
    args=Argument_Parser.parse_args()
    try:
        main()
    except:
        print "Not enough arguments, file not found or an exception occured, try '-h' for help"
        