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


''' Key for itertools groupby
    Alters flag when sequence changes from one condition to another
    Input: Sequence of characters with some alphabet and a trigger condition
    Output: flag: 0 or 1
'''

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

''' kmersWithAmbigIndex
Return list of kmers, indices of k-mers without ambiguity, and indices of those
with ambiguity
Input: string
Output: List of string, list of indices, list of indices
'''
def kmersWithAmbigIndex(tSeq,k):
    O = [tSeq[i:i+k].upper() for i in range(0,len(tSeq)-k+1)]
    O = [o for o in O if 'N' not in o]
    oIdx = [i for i in range(0,len(tSeq)-k+1) if 'N' not in tSeq[i:i+k]]
    nBP = [i for i in range(0,len(tSeq)-k+1) if 'N' in tSeq[i:i+k]]
    nBP = list(zip(nBP,['N']*len(nBP)))
    return O, oIdx, nBP

''' LLR
Return log-likelihood ratio between two models in HMM for + k-mers
Input: sequnce of hits, value of k, k-mer frequencies in HMM emmission matrix
Output: Array of scores for each hit
'''

def LLR(hits,k,E):
    arr = np.zeros(len(hits))
    for i,hit in enumerate(hits):
        scorePos,scoreNeg=0,0
        for j in range(len(hit)-k+1):
            kmer=hit[j:j+k]
            scorePos += E['+'][kmer]
            scoreNeg += E['-'][kmer]
        llr = scorePos-scoreNeg
        arr[i] = llr
    return arr

''' calculateSimilarity
Run several functions including viterbi algorithm, log-likelihood, and generate output dataframes
Input: fasta file information
Output: dataframe object
'''
def calculateSimilarity(data):
    tHead,tSeq = data
    O,oIdx,nBP = kmersWithAmbigIndex(tSeq,k)
    A,E,states,pi= hmm['A'],hmm['E'],hmm['states'],hmm['pi']
    # Viterbi algorithm
    bTrack = corefunctions.viterbi(O,A,E,states,pi)
    #Zip the indices of unambig k-mers with their HMM state labels
    coordBTrack = list(zip(oIdx,bTrack)) # [(1,'-'),(2,'+',...(n,'+'))]
    mergedTrack = coordBTrack + nBP # add back in ambig locations
    mergedTrack.sort(key=itemgetter(0)) # sort master list by index
    hmmTrack = [i[1] for i in mergedTrack] # fetch just state label from mergedTrack ['-','+',...,'+']
    groupedHits = groupHMM(hmmTrack) # ['-----','++++++++++','-','++++','------------']
    idx = 0
    indexGroupHits = []
    '''
    Loop below formats the hmm output as such:

    [([0,1,2]),'---'),([3,4],'++'),([5],'-'),...]

    Grouping HMM states with their correct index in the list of k-mers

    '''
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
            seqHits.append(tSeq[start:end])
    if seqHits:
        starts = [c.split(':')[0] for c in seqHitCoords]
        ends = [c.split(':')[1] for c in seqHitCoords]
        info = list(zip(seqHits,starts,ends))
        dataDict = dict(zip(list(range(len(seqHits))),info))
        df = pd.DataFrame.from_dict(dataDict,orient='index')

        #calculate log-likelihood ratio of k-mers in the + model vs - model
        df['Score'] = LLR(seqHits,k,E)

        df.columns = ['Sequence','Start','End','Score']
        df.sort_values(by='Score',inplace=True,ascending=False)
        df.reset_index(inplace=True)
        fa = df['Sequence']
        df = df[['Start','End','Score']]
        return tHead,[df,fa]
    else:
        return tHead,None


parser = argparse.ArgumentParser()
parser.add_argument("--model",type=str,help='Path to directory containing hmm models from train.py')
parser.add_argument('--db',type=str,help='Path to fasta file with sequences to calculate similarity score')
parser.add_argument('--prefix',type=str,help='String, Output file prefix;default=None')
#parser.add_argument('--bkg',type=str,help='Path to fasta file from which to calculate background nucleotide frequencies, if not passed default is uniform',default=None)
parser.add_argument('-a',type=str,help='String, Alphabet to generate k-mers (e.g. ATCG); default=ATCG',default='ATCG')
parser.add_argument('-n',type=int,help='Integer 1 <= n <= max(cores), Number of processor cores to use; default = 1. This scales with the number of sequence comparisons in --db',default=1)
parser.add_argument('--fasta',action='store_true')

args = parser.parse_args()
alphabet = [letter for letter in args.a]
outLog = open('./log.txt','w')

if os.path.isdir(args.model):
    models = [f+'/' for f in glob.iglob(f'{args.model}*')]
else:
    models = [args.model]

# if args.bkg:
#     bkgFa = Reader(args.bkg)
#     bkgSeqs = bkgFa.get_seqs()
#     probMap = corefunctions.nucContent(bkgSeqs,args.a)
#     outLog.write(f'Background Frequencies: {probMap}')
# elif not args.bkg:
#     probMap = {'A':.25,'T':.25,'C':.25,'G':.25}

for model in models:
    kVals = [f[-1] for f in glob.iglob(f'{model}*')]
    for kVal in kVals:
        kDir = model+kVal+'/'
        modelName = model.split('/')[-2]
        hmm = pickle.load(open(kDir+'hmm.mkv','rb'))
        # Explicitly determine k from the size of the log matrix and the size of the alphabet used to generate it
        k = int(log(len(hmm['E']['+'].keys()),len(args.a)))
        kmers = [''.join(p) for p in itertools.product(alphabet,repeat=k)]
        outLog.write('\nGenerating model of score distribution')
        target = Reader(args.db)
        targetSeqs,targetHeaders = target.get_seqs(),target.get_headers()
        targetMap = defaultdict(list)
        outLog.write('\nScanning database sequences')
        with pool.Pool(args.n) as multiN:
            jobs = multiN.starmap(calculateSimilarity,product(*[list(zip(targetHeaders,targetSeqs))]))
            dataDict = dict(jobs)
        outLog.write('\nDone')
        with open(f'./{args.prefix}_{modelName}_{k}.txt','w') as outfile:
            for h,df in dataDict.items():
                if df:
                    outfile.write(f'{h}\n')
                    outfile.write(df[0].to_string())
                    outfile.write(f'\n')

        if args.fasta:
            with open(f'./{args.prefix}_{modelName}_{k}.fa','w') as outfasta:
                for h,df in dataDict.items():
                    if df:
                        for i,seq in enumerate(df[1]):
                            outfasta.write(f'{h}_hit{i}\n')
                            outfasta.write(f'{seq}\n')

outLog.close()
