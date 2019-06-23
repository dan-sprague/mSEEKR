import numpy as np
import itertools
from itertools import product
from seekr.kmer_counts import BasicCounter
from seekr.fasta_reader import Reader
from tqdm import tqdm as tqdm
from seekr.fasta_reader import Reader
from itertools import groupby
from collections import defaultdict
import seaborn as sns
import matplotlib.pyplot as plt

def score(seq, k, likelihood,alphabet):
    tot=0
    kmers = [seq[i:i+k] for i in range(len(seq)-k+1)]
    nextState = dict(zip(alphabet,range(len(alphabet))))
    currState = dict(zip([''.join(p) for p in product(alphabet,repeat=k-1)],range(4**(k-1))))
    for kmer in kmers:
        if ('N' not in kmer) and ('$' not in kmer):
            i, j = currState[kmer[:k-1]], nextState[kmer[-1]]
            tot += likelihood[i, j]
    return tot

def transitionMatrix(kmers,k,alphabet):
    states = np.zeros((4**(int(k)-1), 4), dtype=np.float64)
    stateKmers = [''.join(p) for p in itertools.product(alphabet,repeat=k-1)]
    for i, currState in enumerate(stateKmers):
        tot = 0
        for j, nextState in enumerate(alphabet):
            count = kmers[currState+nextState]
            tot += count
        if tot > 0:
            for j, nextState in enumerate(alphabet):
                states[i, j] = kmers[currState+nextState] / float(tot)
    return states

def nucContent(nullSeqs,alphabet):
    countDict = dict(zip(alphabet,np.zeros(len(alphabet))))
    nullSeqs = ''.join(nullSeqs)
    for nt in nullSeqs:
        if nt in countDict.keys():
            countDict[nt]+=1

def trainModel(fasta,k,kmers,alphabet):
    q = BasicCounter(fasta,k=k,mean=False,std=False,log2=False,alphabet=alphabet)
    q.get_counts()
    #Reverse length normalization of seekr
    qUnNormCounts = q.counts.T*[len(s) for s in q.seqs]/1000
    qCounts = np.rint(qUnNormCounts.T)
    # qCounts+=1
    qCounts = np.mean(qCounts,axis=0)
    currKmers = dict(zip(kmers,qCounts))
    qTransMat = transitionMatrix(currKmers,k,alphabet)
    return qTransMat

def logLTbl(q,null):
    return np.log2(q) - np.log2(null)


def plotTiles(arr,outname,S):
    sns.set_context('talk')
    plt.figure(figsize=(10,6))
    plt.plot(arr)
    plt.axhline(y=S,linestyle='--',color='black')
    plt.tight_layout()
    plt.savefig(outname,bbox_inches='tight')
    plt.clf()
