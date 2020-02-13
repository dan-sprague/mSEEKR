import pandas as pd
import numpy as np
from seekr.fasta_reader import Reader
from itertools import product
from collections import defaultdict
from scipy.stats import pearsonr
from seekr.kmer_counts import BasicCounter as Count
import glob
from os.path import basename
import pickle
import argparse
import multiprocessing
'''
# This version of SEEKR does not calculate the full square correlation matrix (all pairwise) ,
  but rather a 'rectangular' correlation matrix, representing query-to-target comparisons only
  This leaves out target-target correlations and query-to-query correlations
  If the number of targets is quite large this will obviously save a huge amount of time
'''
###########################################################################

###########################################################################


def rectCorr(S,Q):

    S = S-np.mean(S,axis=1,keepdims=True)
    Q = Q-np.mean(Q,axis=1,keepdims=True)
    covariance = S.dot(Q.T)
    sN = np.linalg.norm(S,axis=1)
    qN = np.linalg.norm(Q,axis=1)
    norm = np.outer(sN,qN)

    CORR = covariance/norm
    return CORR.ravel()

def weight(seqLens,CORR):
    totLen = sum(seqLens)
    weights = [i/totLen for i in seqLens]
    weightedSEEKR = CORR.dot(weights)
    return weightedSEEKR

###########################################################################

###########################################################################
parser = argparse.ArgumentParser()
parser.add_argument('--query',type=str,default=None)
parser.add_argument('--parse',type=str,default=None)
parser.add_argument('--ref',type=str,default=None)
parser.add_argument('-k', type=int,)
parser.add_argument('-n', type=int, help='Number of processors,default = number cpus avail - 1',
                    default=multiprocessing.cpu_count()-1)
parser.add_argument('--out',type=str,)
args = parser.parse_args()

kmers = [''.join(p) for p in product('AGTC',repeat=args.k)]
kmer_map = dict(zip(kmers,range(0,4**args.k)))
###########################################################################

###########################################################################
df = pd.read_csv(args.parse,sep='\t',index_col='seqName')
df = df.drop('Unnamed: 0',axis=1)

ref = Count(infasta=args.ref,
           mean=False,std=False,log2=False,k=args.k)
ref.get_counts()

R = ref.counts
R+=1
R = np.log2(R)

mean,std = np.mean(R,axis=0),np.std(R,axis=0)

query = Count(infasta=args.query,
             mean=False,std=False,log2=False,k=args.k)

query.get_counts()
Q = query.counts
Q+=1
Q=np.log2(Q)
Q = (Q - mean)/std
###########################################################################

'''
Calculate SEEKR correlation of the HMM parses to the query
'''

###########################################################################

seqs = df.Sequence.values

parseCounter = Count(k=args.k,mean=False,std=False,log2=False)
parseCounter.seqs = seqs
parseCounter.get_counts()

S = parseCounter.counts+1
S = np.log2(S)
S = (S-mean)/std

CORR = rectCorr(S,Q)

df['SEEKR'] = CORR
df.drop('Sequence',axis=1,inplace=True)
df.to_csv(f'./{args.out}_{args.k}.txt',sep='\t',header=False)
