import numpy as np
import itertools
from itertools import product
from seekr.kmer_counts import BasicCounter
from seekr.fasta_reader import Reader
from tqdm import tqdm as tqdm
from seekr.fasta_reader import Reader
from itertools import groupby
from collections import defaultdict

def classify(seq, k, lrTab,alphabet):
    seq = seq.upper()
    bits = 0
    nucmap = dict(zip(alphabet,range(len(alphabet))))
    rowmap = dict(zip([''.join(p) for p in product(alphabet,repeat=k-1)],range(4**(k-1))))
    for kmer in [seq[i:i+k] for i in range(len(seq)-k+1)]:
        if ('N' not in kmer) and ('$' not in kmer):
            i, j = rowmap[kmer[:k-1]], nucmap[kmer[-1]]
            #print(f'P({kmer[-1]}|{kmer[:-1]})',lrTab[i,j])
            bits += lrTab[i, j]
    return bits

def markov_chain(kmers,k,alphabet):

    conds = np.zeros((4**(int(k)-1), 4), dtype=np.float64)

    margs = np.zeros(4, dtype=np.float64)

    for i, ci in enumerate([''.join(p) for p in itertools.product(alphabet,repeat=k-1)]):

        tot = 0

        for j, cj in enumerate(alphabet):

            count = kmers[ci+cj]
            tot += count

        if tot > 0:

            for j, cj in enumerate(alphabet):

                conds[i, j] = kmers[ci+cj] / float(tot)

    return conds


def trainModel(fasta,k,kmers,alphabet):
    q = BasicCounter(fasta,k=k,mean=False,std=False,log2=False,alphabet=alphabet)
    q.get_counts()
    #Reverse length normalization of seekr
    qUnNormCounts = q.counts.T*[len(s) for s in q.seqs]/1000
    qCounts = np.rint(qUnNormCounts.T)
    qCounts+=1
    qCounts = np.mean(qCounts,axis=0)
    currKmers = dict(zip(kmers,qCounts))
    qTransMat = markov_chain(currKmers,k,alphabet)

    return qTransMat

def logLTbl(q,null):
    return np.log2(q) - np.log2(null)
