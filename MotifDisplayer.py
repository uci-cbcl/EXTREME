#!/usr/bin/env python
from argparse import ArgumentParser
from pylab import *
from numpy import *

"""
The main executable function
"""
def main():
    usage = "usage: %prog <Data file> <Data folder>"
    description = "The program displays the discovered motifs in an Online EM run"
    parser = ArgumentParser(description=description)
    parser.add_argument('datafile', metavar='d', help='Data file')
    parser.add_argument('folder', metavar='f', help='Folder containing motifs')
    args = parser.parse_args()
    data = open(args.datafile, 'r')
    motifs = list()
    for line in data:
        if "Motif" in line and "value" in line and "e+308" not in line:
            motifs.append(int(line.split()[1]))
    data.close()
    n = len(motifs)
    s = int(ceil(sqrt(n)))#I want a square
    print motifs
    for i in range(1,n+1):
        subplot(s,s,i)
        img = imread(args.folder + '/Motif_' + str(motifs[i-1]) + '.png')
        imshow(img)
    show()

if __name__=='__main__':
    main()