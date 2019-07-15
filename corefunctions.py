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
    obs = [seq[i:i+k] for i in range(len(seq)-k+1)]
    stateKmers = [''.join(p) for p in product(alphabet,repeat=k-1)]
    nextState = dict(zip(alphabet,range(len(alphabet))))
    currState = dict(zip(stateKmers,range(4**(k-1))))
    for kmer in obs:
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
    seqs = ''.join(nullSeqs)
    seqs = seqs.upper()
    freq = [seqs.count(nt)/len(seqs) for nt in alphabet]
    return dict(zip(alphabet,freq))


def trainModel(fasta,k,kmers,alphabet):
    q = BasicCounter(fasta,k=k,mean=False,std=False,log2=False,alphabet=alphabet)
    q.get_counts()
    #Reverse length normalization of seekr
    qUnNormCounts = q.counts.T*[len(s) for s in q.seqs]/1000
    qCounts = np.rint(qUnNormCounts.T)
    qCounts+=1
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

def HMM(S,k,alphabet):
    model,null = S[0],S[1]
    hmmDict = {}
    kmers = [''.join(p) for p in product(alphabet,repeat=k)]
    q = BasicCounter(model,k=k,mean=False,std=False,log2=False,alphabet=alphabet)
    q.get_counts()
    #Reverse length normalization of seekr
    qUnNormCounts = q.counts.T*[len(s) for s in q.seqs]/1000
    qCounts = np.rint(qUnNormCounts.T)
    qCounts+=1
    qCounts = np.mean(qCounts,axis=0)
    hmmDict['model'] = np.log2(qCounts/np.sum(qCounts))
    q = BasicCounter(null,k=k,mean=False,std=False,log2=False,alphabet=alphabet)
    q.get_counts()
    #Reverse length normalization of seekr
    qUnNormCounts = q.counts.T*[len(s) for s in q.seqs]/1000
    qCounts = np.rint(qUnNormCounts.T)
    qCounts+=1
    qCounts = np.mean(qCounts,axis=0)
    #hmmDict['null'] = np.log2(qCounts/sum(qCounts))
    hmmDict['null'] = np.log2(qCounts/np.sum(qCounts))
    states = ('Model','Null')
    pi = {'Model':np.log2(.5),'Null':np.log2(.5)}
    A = {'Model':{'Model':np.log2(.99),'Null':np.log2(.01)},'Null':{'Model':np.log2(.001),'Null':np.log2(.999)}}
    E = {'Model': dict(zip(kmers,hmmDict['model'])),'Null':dict(zip(kmers,hmmDict['null']))}
    return A,E,states,pi

def viterbi(O,A,E,states,pi):

    uk=[{}]
    ukprev = [{}]
    N = len(O)
    for state in states:
        uk[0][state]=pi[state]+E[state][O[0]]
        ukprev[0][state] = None
    for n in range(1,N):
        uk.append({})
        ukprev.append({})
        for state in states:
            prevSelState = states[0]
            currMaxProb = A[state][prevSelState] + uk[n-1][prevSelState]
            for pState in states[1:]:
                currProb = A[state][pState] + uk[n-1][pState]
                if currProb > currMaxProb:
                    currMaxProb = currProb
                    prevSelState = pState
            max_prob = currMaxProb + E[state][O[n]]
            uk[n][state] = max_prob
            ukprev[n][state] = prevSelState

    z = max(uk[-1],key=uk[-n].get)
    prev = ukprev[-1][z]
    backtrack = [z]
    for n in range(N-2,-1,-1):
        backtrack.append(ukprev[n+1][prev])
        prev = ukprev[n+1][prev]
    backtrack = backtrack[::-1]
    return backtrack
