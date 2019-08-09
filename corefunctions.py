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
    LLR=0
    obs = [seq[i:i+k] for i in range(len(seq)-k+1)]
    stateKmers = [''.join(p) for p in product(alphabet,repeat=k-1)]
    nextState = dict(zip(alphabet,range(len(alphabet))))
    currState = dict(zip(stateKmers,range(4**(k-1))))
    for kmer in obs:
        if ('N' not in kmer) and ('$' not in kmer):
            i, j = currState[kmer[:k-1]], nextState[kmer[-1]]
            LLR += likelihood[i, j]
    return LLR

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

def HMM(qCounts,nCounts,k,alphabet,m,n):
    kmers = [''.join(p) for p in itertools.product(alphabet,repeat=k)]
    hmmDict = {}
    countArr = np.array(list(qCounts.values()))
    hmmDict['+'] = np.log2(countArr/np.sum(countArr))

    countArr = np.array(list(nCounts.values()))
    hmmDict['-'] = np.log2(countArr/np.sum(countArr))
    states = ('+','-')
    pi = {'+':np.log2(.5),'-':np.log2(.5)}
    A = {'+':{'+':np.log2(m),'-':np.log2(1-m)},'-':{'+':np.log2(1-n),'-':np.log2(n)}}
    E = {'+': dict(zip(kmers,hmmDict['+'])),'-':dict(zip(kmers,hmmDict['-']))}
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

def baumWelch(O,A,pi,states,E):
    ai = [{}]
    N = len(O)
    for state in states:
        ai[0][state] = pi[state]+E[state][O[0]]
    for n in range(1,N):
        ai.append({})
        for state in states:
            ai[n][state] = E[state][O[n]]
            s = 0
            for j in states:
                s+=ai[n-1][j]+A[j][state]
            ai[n][state]+=s
    BT = [{}]
    for state in states:
        BT[0][state] = np.log2(1)
    backO = O[::-1]
    for n in range(N-1,-1,-1):
        BT.append({})
        s = 0
        for state in states:
            for nextState in states:
                s+=BT[n-1][nextState] + A[state][nextState]+E[state][backO[n-1]]
            BT[n][state] = s
    return ai
