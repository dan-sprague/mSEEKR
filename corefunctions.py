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
from scipy.special import logsumexp
'''

SCORE and transitionMatrix are OLD markov chain functions
Keeping here in case it proves useful

'''

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

# calculate nucleotide content of a sequence... unused now
def nucContent(nullSeqs,alphabet):
    seqs = ''.join(nullSeqs)
    seqs = seqs.upper()
    freq = [seqs.count(nt)/len(seqs) for nt in alphabet]
    return dict(zip(alphabet,freq))

# Calculate the log2 odds table between two matrices
def logLTbl(q,null):
    return np.log2(q) - np.log2(null)


# Plot tiled markov chain plots... unused now
def plotTiles(arr,outname,S):
    sns.set_context('talk')
    plt.figure(figsize=(10,6))
    plt.plot(arr)
    plt.axhline(y=S,linestyle='--',color='black')
    plt.tight_layout()
    plt.savefig(outname,bbox_inches='tight')
    plt.clf()

'''
HMM: Generate dictionaries containing information necessary to perform the viterbi algorithm

Inputs: qCounts - query count dictionary
        nCounts - null count dictionary
        k - value of k
        alphabet - alphabet (ATCG)
        m: + to + transition probability
        n: - to - transition probability
Returns:    A - Dictionary, Hidden state transition matrix
            E - Dictionary, Emission matrix, kmer counts
            state - tuple, list of states (+,-)
            pi - Dictionary, initial probability of being in + or -


'''

def HMM(qCounts,nCounts,k,alphabet,m,n):
    kmers = [''.join(p) for p in itertools.product(alphabet,repeat=k)]
    hmmDict = {}
    countArr = np.array(list(qCounts.values()))
    # Convert raw counts to frequencies, then log probability
    hmmDict['+'] = np.log2(countArr/np.sum(countArr))

    countArr = np.array(list(nCounts.values()))
    # Convert raw counts to frequencies, then log probability
    hmmDict['-'] = np.log2(countArr/np.sum(countArr))
    states = ('+','-')
    pi = {'+':np.log2(.5),'-':np.log2(.5)}
    A = {'+':{'+':np.log2(m),'-':np.log2(1-m)},'-':{'+':np.log2(1-n),'-':np.log2(n)}}
    E = {'+': dict(zip(kmers,hmmDict['+'])),'-':dict(zip(kmers,hmmDict['-']))}
    return A,E,states,pi


'''
Viterbi: Calculate the most likely sequence of hidden states given observed sequence, transition matrix, and emission matrix

Inputs: O - list, observed sequence of k-mers
        A - dictionary, transition matrices of hidden states
        E - dictionary, emission matrices of hidden states
        states - list, hidden states (+,-)
        m: + to + transition probability
        n: - to - transition probability
        pi - Dictionary, initial probability of being in + or -
Returns:    backTrack - list, sequence of hidden states


'''
def viterbi(O,A,E,states,pi):

    # Initialize list of dictionaries for the current step
    # and ukprev, which tracks the state transitions that maximize the 'viterbi function'
    uk=[{}]
    ukprev = [{}]
    N = len(O)
    # calculate initial probabilities in each state given the first kmer
    for state in states:
        uk[0][state]=pi[state]+E[state][O[0]]
        ukprev[0][state] = None # previous state does not exist, set to None
    # Loop through observed sequence
    # For each state, calculate the cumulative probability recursively
    # Store the state transition that maximizes this probability in ukprev for each current state
    for n in range(1,N):
        uk.append({})
        ukprev.append({})
        for state in states:
            prevSelState = states[0] # this is just an arbitrary choice to start checking at the start of the list
            currMaxProb = A[state][prevSelState] + uk[n-1][prevSelState] # probability function
            for pState in states[1:]: # now check the other states...
                currProb = A[state][pState] + uk[n-1][pState]
                if currProb > currMaxProb: # if larger then the first one we checked, we have a new winner, store and continue loop and repeat
                    currMaxProb = currProb
                    prevSelState = pState
            # The emission probability is constant so add at the end rather than in the loop
            max_prob = currMaxProb + E[state][O[n]]
            # save the cumalitive probability for each state
            uk[n][state] = max_prob
            # save the previous state that maximized this probability above
            ukprev[n][state] = prevSelState

    z = max(uk[-1],key=uk[-n].get) # retrieve the state associated with largest log probability
    prev = ukprev[-1][z] # get the state BEFORE "z" above that maximized z
    backtrack = [z,prev] # start our backtrack with knowns
    # Loop through BACKWARDS, getting the previous state that yielded the 'current' max state
    for n in range(N-2,-1,-1):
        backtrack.append(ukprev[n+1][prev]) # n+1 is the "current" state as we're going backwards, ukprev[n+1][prev] returns the state that maximized
        prev = ukprev[n+1][prev]
    fwdTrack = backtrack[::-1] # reverse the order
    return fwdTrack

'''
incomplete baumwelch implementation

'''
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

def fwd(O,A,pi,states,E,k,alphabet):
    a = [{}]
    N = len(O)
    for state in states:
        a[0][state] = pi[state]+E[state][O[0]]
    for n in range(1,N):
        a.append({})
        for state in states:
            P=[]
            naiveP = []
            for pState in states:
                P.append(a[n-1][pState]+A[state][pState] + E[state][O[n]])
            P = logsumexp(P)
            a[n][state] = P
    fwdP = []
    for state in states:
        fwdP.append(a[-1][state])
    fwdP = logsumexp(fwdP)
    return fwdP
