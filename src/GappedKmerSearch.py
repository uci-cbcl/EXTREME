#!/usr/bin/env python
from argparse import ArgumentParser
import string
import sequence
from math import sqrt
from numpy import *

def count_seqs_with_words(seqs, halflength, ming, maxg):
    gapped_seqs_with_words = {}#each key is the gap length
    for g in range(ming,maxg+1):
        seqs_with_words = {}#the current dictionary for the gaps
        gapped_seqs_with_words[g] = seqs_with_words
        print "Looking for k-mers of gap", g
        w = g+2*halflength
        for seq in seqs:
            slen = len(seq)
            for i in range(0, slen-w+1):
                word = seq[i : i+w]
                # skip word if either half-site contains an ambiguous character
                if "N" in word[0:halflength] or "N" in word[-halflength:]:
                    continue
                #convert word to a gapped word. Only the first and last halflength letters are preserved
                word = word[0:halflength] + g*"N" + word[-halflength:]
                update_seqs_with_words(seqs_with_words, word)
    return gapped_seqs_with_words

def update_seqs_with_words(seqs_with_words, word):
    #use the lower alphabet word for rc
    word = min(word, get_rc(word))
    if seqs_with_words.has_key(word):#word has been encountered before, add 1
        seqs_with_words[word] += 1
    else:#word has not been encountered before, create new key
        seqs_with_words[word] = 1

def get_rc(re):
    """ Return the reverse complement of a DNA RE.
    """
    return re.translate(string.maketrans("ACGTURYKMBVDHSWN", "TGCAAYRMKVBHDSWN"))[::-1]

def get_zscores(pos_seq_counts,neg_seq_counts):
    results = {}
    for g in pos_seq_counts:
        results[g] = {}
        for word in pos_seq_counts[g]:
            p = pos_seq_counts[g][word]
            if (neg_seq_counts[g].has_key(word)):
                n = neg_seq_counts[g][word]
            else:
                n = 1
            zscore = 1.0*(p - n)/sqrt(n)
            results[g][word] = zscore
    return results

#returns the words in order, from largest to smallest, by z-scores
def sorted_zscore_keys(zscores):
    keys = {}
    for g in zscores:
        keys[g] = sorted(zscores[g], key=zscores[g].__getitem__, reverse=True)
    return keys

def find_kmers(pos_seqs, neg_seqs, halflength, ming, maxg, minsites, zthresh, outputfile):
    print 'Counting words in positive sequences...'
    pos_seq_counts = count_seqs_with_words(pos_seqs, halflength, ming, maxg)
    print 'Counting words in negative sequences...' 
    neg_seq_counts = count_seqs_with_words(neg_seqs, halflength, ming, maxg)
    print 'Calculating z-scores...'
    zscores = get_zscores(pos_seq_counts,neg_seq_counts)
    print 'Sorting keys by z-scores...'
    sorted_keys = sorted_zscore_keys(zscores)
    output = open(outputfile,'w')
    """
    for g in range(ming,maxg+1):
        subplot(len(range(ming,maxg+1)),1,g+1)
        hist(zscores[g].values(),100,range=(-50,50))
    show()    
    """
    print 'Writing to results to', outputfile
    for g in range(ming,maxg+1):
        numkmers = len(sorted_keys[g])
        z_std = std(zscores[g].values())
        print "Writing k-mers with gap", g
        print "Gap", g, "has standard deviation z-score of", z_std
        for k in range(numkmers):#Keep iterating thru k-mers until z-scores below threshold
            key = sorted_keys[g][k]
            zscore = zscores[g][key]
            corrected_zscore = zscore/z_std
            if corrected_zscore < zthresh or pos_seq_counts[g][key] < minsites:
                print "Gap", g, "had", k, "k-mers above the threshold that are above the minimum number of sites"
                break
            pos_sites = pos_seq_counts[g][key]
            if neg_seq_counts[g].has_key(key):
                neg_sites = neg_seq_counts[g][key]
            else:
                neg_sites = 0
            output.write(str(key)+"\t"+str(pos_sites)+"\t"+str(neg_sites)+"\t"+str(corrected_zscore)+"\t"+str(zscore)+"\n")
    output.close()

"""
The main executable function
"""
def main():
    usage = "usage: %prog [options] <input FASTA> <negative FASTA>"
    description = "The program performs a DREME-like search for gapped k-mers"
    parser = ArgumentParser(description=description)
    parser.add_argument('fastafile', metavar='f', help='FASTA file containing the sequences')
    parser.add_argument('negativefile', metavar='n', help='FASTA file containing the negative sequences')
    parser.add_argument('outputfile', metavar='o', help='Output file')
    parser.add_argument("-w", "--width", dest="width", help="Width of the motif to search for. This makes the program only search for a motif of this width. Beware if greater than 8", type=int, default=0)
    parser.add_argument("-ming", dest="mingap", help="Minimum gap of k-mer to search for. Default: 0", type=int, default=0)
    parser.add_argument("-maxg", dest="maxgap", help="Maximum gap of k-mer to search for. Default: 12", type=int, default=10)
    parser.add_argument("-l", dest="halflength", help="Number of non-degenerate letters per half-site. Total number of non-degenerate letters is twice this number. Default: 4", type=int, default=4)
    parser.add_argument("-minw", dest="minwidth", help="Minimum width of the motif to search for. The default is 3, which is the width of the smallest core motif.", type=int, default=3)
    parser.add_argument("-maxw", dest="maxwidth", help="Maximum width of the motif to search for. This program does one refinement at this width (if greater than 8), and then picks the most significant short-mer. Default: 8", type=int, default=8)
    parser.add_argument("-mink", dest="mink", help="Minimum width of the core to search for. The default is 3, which is the width of the smallest core motif.", type=int, default=3)
    parser.add_argument("-maxk", dest="maxk", help="Maximum width of the core to search for. Default: 8", type=int, default=8)
    parser.add_argument("-z", "--zthresh", dest="zthresh", help="Corrected z-score threshold. Default: 5", type=float, default=5)
    parser.add_argument("-minsites", "--minsites", dest="minsites", help="Minimum number of sites for a k-mer to be included. Default: 10", type=int, default=10)    
    args = parser.parse_args()
    pos_seq_file_name = args.fastafile
    neg_seq_file_name = args.negativefile
    print 'Reading positive sequence file...'    
    pos_seqs = sequence.convert_ambigs(sequence.readFASTA(pos_seq_file_name, None, True))
    print 'Reading negative sequence file...'
    neg_seqs = sequence.convert_ambigs(sequence.readFASTA(neg_seq_file_name, None, True))
    halflength = args.halflength
    ming = args.mingap
    maxg = args.maxgap
    zthresh = args.zthresh
    minsites = args.minsites
    find_kmers(pos_seqs, neg_seqs, halflength, ming, maxg, minsites, zthresh, args.outputfile)

if __name__=='__main__':
    main()
