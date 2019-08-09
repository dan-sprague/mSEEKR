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
import pickle

parser = argparse.ArgumentParser()
parser.add_argument("--query",type=str,help='Path to kmer count file for sequences of interest (e.g. functional regions of a ncRNA)')
parser.add_argument('--null', type=str,help='Path to kmer count file that compose null model (e.g. transcriptome, genome, etc.)')
parser.add_argument('--qT',type=float,help='Probability of query to query transition',default=.999)
parser.add_argument('--nT',type=float,help='Probability of null to null transition',default=.9999)
parser.add_argument('--qPrefix',type=str,help='String, Output file prefix;default=None',default='query')
parser.add_argument('--nullPrefix',type=str,help='String, Output file prefix;default=None',default='null')
parser.add_argument('--dir',type=str,help='Output directory',default='./markovModels/')
parser.add_argument('-k',type=str,help='Comma delimited string of possible k-mer values,must be found in the k-mer count file',default='2,3,4')
parser.add_argument('-a',type=str,help='String, Alphabet to generate k-mers (e.g. ATCG); default=ATCG',default='ATCG')
args = parser.parse_args()
kVals = [int(i) for i in args.k.split(',')]

newDir = f'{args.dir}{args.qPrefix}_{args.nullPrefix}/'
if not os.path.exists(newDir):
    os.mkdir(newDir)
S = [args.query,args.null]
alphabet = [letter for letter in args.a]

qCount = pickle.load(open(args.query,'rb'))
nCount = pickle.load(open(args.null,'rb'))

for k in kVals:
    if (k in qCount.keys()) and (k in nCount.keys()):
        qKCount = qCount[k]
        nKCount = nCount[k]
        kDir = newDir+f'{k}/'
        if not os.path.exists(kDir):
            os.mkdir(kDir)
        A,E,states,pi = corefunctions.HMM(qKCount,nKCount,k,args.a,args.qT,args.nT)
        kmers = [''.join(p) for p in itertools.product(alphabet,repeat=k)]
        queryMkv = corefunctions.transitionMatrix(qKCount,k,alphabet)
        nullMkv = corefunctions.transitionMatrix(nKCount,k,alphabet)
        lgTbl = corefunctions.logLTbl(queryMkv,nullMkv)
    else:
        print(f'Missing {k}-mer counts in count file... skipping')

    np.savetxt(f'{kDir}logtbl.mkv',lgTbl)
    pickle.dump({'A':A,'E':E,'pi':pi,'states':states},open(f'{kDir}hmm.mkv','wb'))
