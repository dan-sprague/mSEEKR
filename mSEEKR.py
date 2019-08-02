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
''' Key for itertools groupby
    Alters flag when sequence changes from annotated from unannotated or vice versa
    Input: genomic sequence
    Output: flag: 0 or 1
'''

class Key(object):
    def __init__(self):
        self.is_nt,self.flag,self.prev = ['-'],[0,1],None
    def __call__(self,e):
        ebool = any(x in self.is_nt for x in e)
        if self.prev:
            prevbool = any(x in self.is_nt for x in self.prev)
        else:
            prevbool = None
        if prevbool != ebool:
            self.flag = self.flag[::-1]
        self.prev = e
        return self.flag[0]

''' fastaN
Return a list of lists separating HMM state labels
Input: String
Output: List of lists
Output example: ['---','++','------','+','-------',...]
'''



def groupHMM(seq):
    return [''.join(list(g)) for k,g in itertools.groupby(seq,key=Key())]

def calculateSimilarity(data):
    tHead,tSeq = data
    O = [tSeq[i:i+k].upper() for i in range(0,len(tSeq)-k+1)]
    O = [o for o in O if 'N' not in o]
    A,E,states,pi= hmm['A'],hmm['E'],hmm['states'],hmm['pi']
    bTrack = corefunctions.viterbi(O,A,E,states,pi)
    groupedHits = groupHMM(bTrack)
    idx = 0
    haha = []
    for i,group in enumerate(groupedHits):
        haha.append([])
        for kmer in group:
            haha[i].append(idx)
            idx+=1
    hits = list(zip(haha,groupedHits))
    seqHits = []
    seqHitCoords = []
    for group in hits:
        if '+' in group[1]:
            start,end = group[0][0],group[0][-1]+k+1
            if end-start >= 25:
                seqHitCoords.append(f'{start}:{end}')
                seqHits.append(tSeq[start:end])


    seqScores = np.array([corefunctions.score(tile,k,lgTbl,alphabet) for tile in seqHits])

    dataDict = dict(zip(seqHitCoords,seqHits))
    df = pd.DataFrame.from_dict(dataDict,orient='index')
    df['Score'] = seqScores
    df.columns = ['Sequence','Score']
    df.index.name = 'bp'
    df.sort_values(by='Score',inplace=True,ascending=False)
    #df = df[df['Score']>1]

    return tHead,[df,len(tSeq)]


parser = argparse.ArgumentParser()
parser.add_argument("--model",type=str,help='Path to directory containing .mkv files or path to a single .mkv file;default=./markovModels/',default='./markovModels/')
parser.add_argument('--db',type=str,help='Path to fasta file with sequences to calculate similarity score')
parser.add_argument('--prefix',type=str,help='String, Output file prefix;default=None')
parser.add_argument('--bkg',type=str,help='Path to fasta file from which to calculate background nucleotide frequencies, if not passed default is uniform',default=None)
parser.add_argument('-a',type=str,help='String, Alphabet to generate k-mers (e.g. ATCG); default=ATCG',default='ATCG')
parser.add_argument('-n',type=int,help='Integer 1 <= n <= max(cores), Number of processor cores to use; default = 1. This scales with the number of sequence comparisons in --db',default=1)


args = parser.parse_args()
alphabet = [letter for letter in args.a]
outLog = open('./log.txt','w')

if os.path.isdir(args.model):
    models = [f+'/' for f in glob.iglob(f'{args.model}*')]
else:
    models = [args.model]

if args.bkg:
    bkgFa = Reader(args.bkg)
    bkgSeqs = bkgFa.get_seqs()
    probMap = corefunctions.nucContent(bkgSeqs,args.a)
    outLog.write(f'Background Frequencies: {probMap}')
elif not args.bkg:
    probMap = {'A':.25,'T':.25,'C':.25,'G':.25}

for model in models:
    kVals = [f[-1] for f in glob.iglob(f'{model}*')]
    for kVal in kVals:
        kDir = model+kVal+'/'
        modelName = model.split('/')[-2]
        lgTbl = np.loadtxt(kDir+'logtbl.mkv')
        hmm = pickle.load(open(kDir+'hmm.mkv','rb'))
        # Explicitly determine k from the size of the log matrix and the size of the alphabet used to generate it
        k = int(log(lgTbl.size,len(args.a)))
        kmers = [''.join(p) for p in itertools.product(alphabet,repeat=k)]
        probs = [probMap[letter] for letter in args.a]
        outLog.write('\nGenerating model of score distribution')
        target = Reader(args.db)
        targetSeqs,targetHeaders = target.get_seqs(),target.get_headers()
        targetMap = defaultdict(list)
        outLog.write('\nScanning database sequences')
        with pool.Pool(args.n) as multiN:
            jobs = multiN.starmap(calculateSimilarity,product(*[list(zip(targetHeaders,targetSeqs))]))
            dataDict = dict(jobs)
        outLog.write('\nDone')
        for h,df in dataDict.items():
            df[0].index.name = f'bp'
            df[0].to_csv(f'./{args.prefix}_{h[1:]}_{modelName}_{k}_{df[1]}.txt',sep='\t')
outLog.close()
