import corefunctions
import coreStats
import argparse
import itertools
import numpy as np
from seekr.fasta_reader import Reader
from scipy.stats import norm
from collections import defaultdict
from multiprocessing import pool
from scipy.stats import gaussian_kde
from itertools import product
import sys
import os
import glob
from math import log
import pickle
import itertools
import pandas as pd
from operator import itemgetter


class Key(object):
    def __init__(self):
        self.is_nt,self.flag,self.prev = ['-','N'],[0,1],None
    def __call__(self,e):
        # Initial True/False  if first char in string is + or -
        ebool = any(x in self.is_nt for x in e)
        # If key exists (self.prev is defined), do true/false check
        # else, set value to false
        if self.prev:
            prevbool = any(x in self.is_nt for x in self.prev)
        else:
            prevbool = None
        # if string goes from - to +, or + to -, swap flag
        if prevbool != ebool:
            self.flag = self.flag[::-1]
        # set previous encountered char, for the next interation of this, to the current value
        self.prev = e
        return self.flag[0]

''' groupHMM
Return a list of strings separating HMM state labels
Input: String
Output: List of lists
Output example: ['---','++','------','+','-------',...]
'''
def groupHMM(seq):
    return [''.join(list(g)) for k,g in itertools.groupby(seq,key=Key())]


'''

=====================
Arguments
=====================

'''

parser = argparse.ArgumentParser()
parser.add_argument("--model",type=str,help='Path to directory containing hmm models from train.py')
parser.add_argument("-k",type=int,help='Value of k to use')
parser.add_argument('--rand',type=int,help='length of random DNA',default=10000)
parser.add_argument('--name',type=str,help='String, Output file prefix;default=None',default='stats')
parser.add_argument('-a',type=str,help='String, Alphabet to generate k-mers (e.g. ATCG); default=ATCG',default='ATCG')
parser.add_argument('--bkg',type=str,help='Comma delimited string of nucleotide background frequencies, provided in the same order as the -a argument,default = uniform',default='.25,.25,.25,.25')
parser.add_argument('-n',type=int,help='Integer 1 <= n <= max(cores), Number of processor cores to use; default = 1. This scales with the number of sequence comparisons in --db',default=1)
args = parser.parse_args()


kVal = args.k
bkgFreq = [float(i) for i in args.bkg.split(',')]
model = args.model

if not model.endswith('/'):
    model +='/'
kDir = model+str(kVal)+'/'
modelName = model.split('/')[-2]
# Check if file exists and open if so, else skip this iteration of the loop
try:
    hmm = pickle.load(open(kDir+'hmm.mkv','rb'))
except:
    print(f'Value of k={kVal} not found, exiting...')
    sys.exit()
# Explicitly determine k from the size of the log matrix and the size of the alphabet used to generate it
k = int(log(len(hmm['E']['+'].keys()),len(args.a)))

E = np.array(list(hmm['E']['+'].values())) - np.array(list(hmm['E']['-'].values()))
a = [i.upper() for i in args.a]
randDNA = ''.join(np.random.choice(a,args.rand,p=bkgFreq))

O = [randDNA[i:i+k].upper() for i in range(0,len(randDNA)-k+1)]
O = [o for o in O if 'N' not in o]
A,E,states,pi= hmm['A'],hmm['E'],hmm['states'],hmm['pi']
bTrack = corefunctions.viterbi(O,A,E,states,pi)

groupedHits = groupHMM(bTrack) # ['-----','++++++++++','-','++++','------------']
idx = 0
indexGroupHits = []

# Loop below formats the hmm output as such:
# [([0,1,2]),'---'),([3,4],'++'),([5],'-'),...]
# Grouping HMM states with their correct index in the list of k-mers
for i,group in enumerate(groupedHits):
    indexGroupHits.append([])
    for kmer in group:
        indexGroupHits[i].append(idx)
        idx+=1
hits = list(zip(indexGroupHits,groupedHits))
seqHits = []
seqHitCoords = []
for group in hits:
    if '+' in group[1]:
        start,end = group[0][0],group[0][-1]+k #convert k-mer coord to bp coord
        seqHitCoords.append(f'{start}:{end}')
        seqHits.append(randDNA[start:end])
starts = np.array([int(c.split(':')[0]) for c in seqHitCoords])
ends = np.array([int(c.split(':')[1]) for c in seqHitCoords])
fwdPs = []
for hit in seqHits:
    O = [hit[i:i+k].upper() for i in range(0,len(hit)-k+1)]
    O = [o for o in O if 'N' not in o]
    '''
    forward algorithm to calculate log P(O|HMM)
    '''
    fP = corefunctions.fwd(O,A,pi,states,E,k,args.a)
    fwdPs.append(fP)

fwdPs = np.array(fwdPs)

np.savetxt('./scoreDist.txt',fwdPs)
