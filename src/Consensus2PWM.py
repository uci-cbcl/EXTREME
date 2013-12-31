#!/usr/bin/env python
from argparse import ArgumentParser
from numpy import *
from collections import OrderedDict

def make_PWMs(consensus_filename, output_filename):
    output = open(output_filename,'w')
    #zscores_dictionary = get_zscores_dictionary(counts_filename)
    #contains the lower alphabet k-mer and its z-score
    PWM_dict = get_PWM_dictionary(consensus_filename)
    keys = PWM_dict.keys()
    for key in keys:
        pwm = PWM_dict[key][0]
        consensus_sequence = PWM_dict[key][1]
        output.write('>'+key+"\t"+consensus_sequence+"\n")
        output.write(array_str(pwm,precision=6).replace('[','').replace(']','').replace('0. ', '0.000000').replace('1. ', '1.000000').replace('0.\n', '0.000000\n').replace('1.\n', '1.000000\n').replace('\n  ','\n')[1:]+"\n\n")
    output.close()

def get_PWM_dictionary(consensus_filename):
    PWM_dict = OrderedDict()
    consensusfile = open(consensus_filename,'r')
    #iterate through consensus
    consensus_lines = consensusfile.readlines()
    i = 0
    while i < len(consensus_lines):
        consensus_line = consensus_lines[i]
        if '>' in consensus_line:
            parts = consensus_line.split()
            key = parts[0][1:]
            n = int(parts[1])
            PWM_dict[key] = get_PWM(consensus_lines[i+1:i+1+n])
            i += n + 2
    consensusfile.close()
    return PWM_dict

"""
Converts the lines following a name (eg. >cluster1) into a PWM. Each
count is weighted by the z-score for a k-mer. Function also returns an
improved consensus sequence.
"""
def get_PWM(lines):
    parts_list = [l.split() for l in lines]
    width = max([len(p[0]) for p in parts_list])
    consensus_parts = []
    dna_ambig_dict = {
        'A' : 'A',
        'C' : 'C',
        'G' : 'G',
        'T' : 'T',
        'AG' : 'R',
        'CT' : 'Y',
        'GT' : 'K',
        'AC' : 'M',
        'CG' : 'S',
        'AT' : 'W',
        'CGT' : 'B',
        'AGT' : 'D',
        'ACT' : 'H',
        'ACG' : 'V',
        'ACGT' : 'N'
        }
    for w in range(width):
        consensus_parts.append(set())
    pwm = zeros((width,4))
    #making the pwm
    for p in parts_list:
        #zscore_key = p[1]
        #zscore = zscores_dictionary[zscore_key]
        zscore = float(p[4])
        kmer = p[0] + '.'*(width - len(p[0]))
        #iterate through the kmer letters and add to the pwm according to each letter
        for k in range(width):
            letter = kmer[k]
            if letter =='.' or letter == 'N':
                pwm[k] += zscore*array([0.23,0.27,0.27,0.23])
            elif letter == 'A':
                consensus_parts[k].add('A')
                pwm[k] += zscore*array([0.88,0.04,0.04,0.04])
            elif letter == 'C':
                consensus_parts[k].add('C')
                pwm[k] += zscore*array([0.04,0.88,0.04,0.04])
            elif letter == 'G':
                consensus_parts[k].add('G')
                pwm[k] += zscore*array([0.04,0.04,0.88,0.04])
            elif letter == 'T':
                consensus_parts[k].add('T')
                pwm[k] += zscore*array([0.04,0.04,0.04,0.88])
    #normalize the pwm
    pwm = pwm/pwm.sum(axis=1)[:,newaxis]
    consensus_sequence = ""
    for consensus_part in consensus_parts:
        consensus_part_sorted_list = sorted(consensus_part)
        consensus_part_string = "".join(consensus_part_sorted_list)
        if consensus_part_string == "":
            consensus_part_string = "ACGT"
        consensus_sequence += dna_ambig_dict[consensus_part_string]
    return (pwm, consensus_sequence)

def get_zscores_dictionary(counts_filename):    
    countsfile = open(counts_filename,'r')
    counts_lines = countsfile.readlines()
    zscores_dictionary = {}
    for counts_line in counts_lines:
        parts = counts_line.split()
        zscores_dictionary[parts[0]] = float(parts[-1]) 
    countsfile.close()
    return zscores_dictionary
    
"""
The main executable function
"""
def main():
    usage = "usage: %prog kmercounts consensusfile outputfile"
    description = "The program takes a consensus output to generate a PWM for each consensus group."
    parser = ArgumentParser(description=description)
    parser.add_argument('consensusfile', metavar='c', help='Output from Clustering containing Kmer clusters')
    parser.add_argument('outputfile', metavar='o', help='Output file for results')
    args = parser.parse_args()
    consensus_filename = args.consensusfile
    output_filename = args.outputfile
    make_PWMs(consensus_filename, output_filename)

if __name__=='__main__':
    main()