from seekr.fasta_reader import Reader
import numpy as np
import corefunctions
import argparse
import pickle
from tqdm import tqdm as tqdm
from scipy.special import logsumexp
from collections import defaultdict
import corefunctions
import argparse
from itertools import product
import numpy as np
from seekr.fasta_reader import Reader
from collections import defaultdict
from multiprocessing import pool
import sys
import os
import pickle
from math import log
import pandas as pd
from operator import itemgetter
import time

parser = argparse.ArgumentParser()
parser.add_argument("-k",type=int)
parser.add_argument('--db',type=str,help='Path to fasta file with sequences to calculate similarity score')
parser.add_argument('--prior',type=str)
args = parser.parse_args()


fa = Reader(args.db)
seqss = fa.get_seqs()
h = fa.get_headers()
model = args.prior

for seqs in seqss:
k = args.k
O,oIdx,nBP = corefunctions.kmersWithAmbigIndex(seqs,k)
hmm = pickle.load(open(args.prior,'rb'))
A,E,pi,states = hmm['A'],hmm['E'],hmm['pi'],hmm['states']

for i in tqdm(range(20)):
    a = corefunctions.fwd(O,A,pi,states,E)
    b = corefunctions.bkw(O,A,pi,states,E)
    A = corefunctions.update(a,b,O,states,A,E)
    for i in states:
        marg = logsumexp(list(A[i].values()))
        A[i]['+']-=marg
        A[i]['-']-=marg


bTrack = corefunctions.viterbi(O,A,E,states,pi)
#Zip the indices of unambig k-mers with their viterbi derived HMM state labels
coordBTrack = list(zip(oIdx,bTrack)) # [(1,'-'),(2,'+',...(n,'+'))]
mergedTrack = coordBTrack + nBP # add back in ambig locations
mergedTrack.sort(key=itemgetter(0)) # sort master list by index
hmmTrack = [i[1] for i in mergedTrack] # fetch just state label from mergedTrack ['-','+',...,'+']
groupedHits = corefunctions.groupHMM(hmmTrack) # ['-----','++++++++++','-','++++','------------']

# Return sequences of HMM hits, and their start and end locations in the original sequence
seqHits,starts,ends = corefunctions.formatHits(groupedHits,k,seqs)

df = corefunctions.hitOutput(seqHits,starts,ends,k,E,h,seqs)

df.to_csv('./viterbi_bw.txt',sep='\t',index=False)