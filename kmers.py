import kmers
import argparse
import pickle
import os
from multiprocessing import pool
from itertools import product
from math import ceil
from seekr.fasta_reader import Reader

'''
This program counts k-mers for multiple specified values of k and saves them
to a binary file that countains a dictionary of dictionaries
'''

#Load arguments, see help= for explanation
parser = argparse.ArgumentParser()
parser.add_argument('--fasta', type=str,help='Path to fasta file')
parser.add_argument('--name',type=str,help='name for count file')
parser.add_argument('--dir',type=str,help='directory to save output',default='./')
parser.add_argument('-k',type=str,help='Comma delimited string of possible k-mer values',default='2,3,4')
parser.add_argument('-a',type=str,help='String, Alphabet to generate k-mers (e.g. ATCG); default=ATCG',default='ATCG')
parser.add_argument('-n',type=int,help='Number of CPU cores. Each job corresponds to a value of k, and the program scales well with multiprocessing',default=1)

args = parser.parse_args()

# Read in specified values of k, and the alphabet
kVals = [int(i) for i in args.k.split(',')]
a = args.a.upper()

#SEEKR fasta reader module
F = Reader(args.fasta)
fS = F.get_seqs()

#Join sequences together using $ delimiter character
fString = '$'.join(fS).upper()
lenFString = sum([len(i) for i in fS])

# Need to figure out how to deal with very long fasta files (~ 2-3X the size of the transcriptome in mice)
# if lenFString >= 2147483647:
#     fString='$'.join(fS[::10]).upper()

#Split jobs onto processors and call kmers.pyx cython file
with pool.Pool(args.n) as multiN:
    jobs = multiN.starmap(kmers.main,product(*[[fString],kVals,[a]]))
    dataDict = dict(jobs)

#Save data 
kDir = args.dir
pickle.dump(dataDict,open(f'{kDir}{args.name}.skr','wb'))
