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


''' calculateSimilarity
Run several functions including viterbi algorithm, log-likelihood, and generate output dataframes
Input: fasta file information
Output: dataframe object
'''
def calculateSimilarity(data):
    tHead,tSeq = data
    O,oIdx,nBP = corefunctions.kmersWithAmbigIndex(tSeq,k)
    A,E,states,pi= hmm['A'],hmm['E'],hmm['states'],hmm['pi']
    # Viterbi algorithm
    bTrack = corefunctions.viterbi(O,A,E,states,pi)
    #Zip the indices of unambig k-mers with their HMM state labels
    coordBTrack = list(zip(oIdx,bTrack)) # [(1,'-'),(2,'+',...(n,'+'))]
    mergedTrack = coordBTrack + nBP # add back in ambig locations
    mergedTrack.sort(key=itemgetter(0)) # sort master list by index
    hmmTrack = [i[1] for i in mergedTrack] # fetch just state label from mergedTrack ['-','+',...,'+']
    groupedHits = corefunctions.groupHMM(hmmTrack) # ['-----','++++++++++','-','++++','------------']
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
            seqHits.append(tSeq[start:end])
    starts = np.array([int(c.split(':')[0]) for c in seqHitCoords])
    ends = np.array([int(c.split(':')[1]) for c in seqHitCoords])
    fwdPs = []
    for hit in seqHits:
        O = [hit[i:i+k].upper() for i in range(0,len(hit)-k+1)]
        O = [o for o in O if 'N' not in o]
        '''
        forward algorithm to calculate log P(O|HMM)
        '''
        fP = corefunctions.fwd(O,A,pi,states,E,k,alphabet)
        fwdPs.append(fP) #negative log P val
    # Standard output (hit by hit)
    if (seqHits) and (not args.wt):
        info = list(zip(seqHits,starts,ends))
        dataDict = dict(zip(list(range(len(seqHits))),info))
        df = pd.DataFrame.from_dict(dataDict,orient='index')
        #calculate log-likelihood ratio of k-mers in the + model vs - model
        df['kmerLLR'] = corefunctions.LLR(seqHits,k,E)
        df['fwdLLR'] = fwdPs
        df['seqName'] = tHead
        df.columns = ['Sequence','Start','End','kmerLLR','fwdLLR','seqName']
        df.sort_values(by='fwdLLR',inplace=True,ascending=False)
        df.reset_index(inplace=True)
        fa = df['Sequence']
        df = df[['Start','End','kmerLLR','fwdLLR','seqName','Sequence']]
        return tHead,df

    # Alternative output (transcript by transcript)
    elif (seqHits) and (args.wt):
        sumHits = corefunctions.LLR(seqHits,k,E)
        lens = ends-starts # lengths of individual hits
        df = pd.DataFrame([np.sum(sumHits)]) # sum of hits
        df['totalLenHits'] = (np.sum(lens)) # sum of all hit lengths
        df['fracTranscriptHit'] = df['totalLenHits']/len(tSeq) # fraction of transcript that is hit
        df['longestHit'] = np.max(lens) # longest HMM hit
        df['seqName'] = tHead
        df['sumFwdAlgLogP'] = (np.sum(fwdPs))
        df.columns = ['sumLLR','totalLenHits','fracTranscriptHit','longestHit','seqName','sumFwdLLR']
        return tHead,df
    else:
        return tHead,None


parser = argparse.ArgumentParser()
parser.add_argument("--model",type=str,help='Path to directory containing hmm models from train.py')
parser.add_argument("-k",type=str,help='Comma delimited string, values of k to use')
parser.add_argument('--db',type=str,help='Path to fasta file with sequences to calculate similarity score')
parser.add_argument('--prefix',type=str,help='String, Output file prefix;default=None')
#parser.add_argument('--bkg',type=str,help='Path to fasta file from which to calculate background nucleotide frequencies, if not passed default is uniform',default=None)
parser.add_argument('-a',type=str,help='String, Alphabet to generate k-mers (e.g. ATCG); default=ATCG',default='ATCG')
parser.add_argument('-n',type=int,help='Integer 1 <= n <= max(cores), Number of processor cores to use; default = 1. This scales with the number of sequence comparisons in --db',default=1)
parser.add_argument('--fasta',action='store_true',help='FLAG: print sequence of hit, ignored if --wt is passed')
parser.add_argument('--wt',action='store_true',help='FLAG: If passed, return total log-likelihood over all hits in a fasta entry')
args = parser.parse_args()
alphabet = [letter for letter in args.a]
#Loop over values of k
kVals = args.k.split(',')
args.a = args.a.upper()
model = args.model

if not model.endswith('/'):
    model +='/'
for kVal in kVals:
    kDir = model+kVal+'/'
    modelName = model.split('/')[-2]
    # Check if file exists and open if so, else skip this iteration of the loop
    try:
        hmm = pickle.load(open(kDir+'hmm.mkv','rb'))
    except:
        print(f'Value of k={kVal} not found, skipping...')
        continue
    # Explicitly determine k from the size of the log matrix and the size of the alphabet used to generate it
    k = int(log(len(hmm['E']['+'].keys()),len(args.a)))
    kmers = [''.join(p) for p in itertools.product(alphabet,repeat=k)] # generate k-mers
    target = Reader(args.db)
    targetSeqs,targetHeaders = target.get_seqs(),target.get_headers()
    targetMap = defaultdict(list)
    #Pool processes onto number of CPU cores specified by the user
    with pool.Pool(args.n) as multiN:
        jobs = multiN.starmap(calculateSimilarity,product(*[list(zip(targetHeaders,targetSeqs))]))
        dataDict = dict(jobs)
    #Check if no hits were found
    # if not all(v == None for v in dataDict.values()):
    if not args.wt:
        dataFrames = pd.concat([df for df in dataDict.values() if not None])
        dataFrames['Length'] = dataFrames['End'] - dataFrames['Start']
        dataFrames = dataFrames[['Start','End','Length','fwdLLR','kmerLLR','seqName','Sequence']]
        if not args.fasta:
            dataFrames = dataFrames[['Start','End','Length','fwdLLR','kmerLLR','seqName']]
        dataFrames.sort_values(by='fwdLLR',ascending=False,inplace=True)
        dataFrames.reset_index(inplace=True,drop=True)
        dataFrames.to_csv(f'./{args.prefix}_{modelName}_{k}.txt',sep='\t')
    elif args.wt:
        dataFrames = pd.concat([df for df in dataDict.values() if not None])
        dataFrames = dataFrames[['seqName','sumFwdLLR','sumLLR','totalLenHits','fracTranscriptHit','longestHit']]
        dataFrames.sort_values(by='sumFwdLLR',ascending=False,inplace=True)
        dataFrames.reset_index(inplace=True,drop=True)
        dataFrames.to_csv(f'./{args.prefix}_{modelName}_{k}.txt',sep='\t')
