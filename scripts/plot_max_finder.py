#!/usr/bin/python2.7
'''
Created on 3/06/2013

@author: suu13
'''


import glob
import argparse

def main():
    plot_stranded=[]
    for p_files in glob.glob(str("*"+args.plotextension)):
        print "Plot file name :"+p_files
        with open(p_files) as plot_file:
            plot_lines=plot_file.readlines()
            if len(plot_stranded)==0:
                for _ in range(0,len(plot_lines)):
                    plot_stranded.append([0,0])
            for i in range(0,len(plot_lines)):
                if int(plot_lines[i].split()[0])>plot_stranded[i][0]:
                    plot_stranded[i][0]=int(plot_lines[i].split()[0])
                if int(plot_lines[i].split()[1])>plot_stranded[i][1]:
                    plot_stranded[i][1]=int(plot_lines[i].split()[1])
    with open(args.output,"w") as output_file:
        for i in range(0,len(plot_stranded)):
            output_file.write(str(plot_stranded[i][0])+'\t'+str(plot_stranded[i][1])+'\n')        
            
            
    


if __name__ == '__main__':
    Argument_Parser=argparse.ArgumentParser(prog="plot_max_finder.py")
    Argument_Parser.add_argument('-plotextension',type=str,help="Plot files extension to find max value",required=True)
    Argument_Parser.add_argument('-output',type=str,help="Output file name",required=True)
    args=Argument_Parser.parse_args()
    main()
    pass