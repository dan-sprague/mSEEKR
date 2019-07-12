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

def HMM(S,k,alphabet,O):
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
    obs = [O[i:i+k] for i in range(0,len(O)-k+1)]
    states = ('Model','Null')
    start_p = {'Model':np.log2(.5),'Null':np.log2(.5)}
    trans_p = {'Model':{'Model':np.log2(.8),'Null':np.log2(.2)},'Null':{'Model':np.log2(.1),'Null':np.log2(.9)}}
    emit_p = {'Model': dict(zip(kmers,hmmDict['model'])),'Null':dict(zip(kmers,hmmDict['null']))}
    return obs,states,start_p,trans_p,emit_p

def viterbi(obs,states,start_p,trans_p,emit_p):
    V = [{}]
    for st in states:
        V[0][st] = {'prob':start_p[st]+emit_p[st][obs[0]],'prev':None}
    for t in range(1,len(obs)):
        V.append({})
        for st in states:
            max_tr_prob = V[t-1][states[0]]['prob']+trans_p[states[0]][st]
            prev_st_selected = states[0]
            for prev_st in states[1:]:
                tr_prob = V[t-1][prev_st]['prob']+trans_p[prev_st][st]
                if tr_prob > max_tr_prob:
                    max_tr_prob = tr_prob
                    prev_st_selected = prev_st
            max_prob = max_tr_prob + emit_p[st][obs[t]]
            V[t][st] = {'prob':max_prob,'prev':prev_st_selected}

    opt = []
    max_prob = max(value['prob'] for value in V[-1].values())
    previous = None
    for st,data in V[-1].items():
        if data['prob'] == max_prob:
            opt.append(st)
            previous = st
            break

    for t in range(len(V)-2,-1,-1):
        opt.insert(0,V[t+1][previous]['prev'])
        previous = V[t+1][previous]['prev']
    return opt
    # uk=[{}]
    # ukprev = [{}]
    # N = len(obs)
    # print(N)
    # for state in states:
    #     uk[0][state]=start_p[state]+emit_p[state][obs[0]]
    #     ukprev[0][state] = None
    # for n in range(1,N):
    #     uk.append({})
    #     ukprev.append({})
    #     for state in states:
    #         prevSelState = states[0]
    #         currMaxProb = emit_p[state][obs[n]]+trans_p[state][prevSelState] + uk[n-1][state]
    #         for pState in states[1:]:
    #             currProb = emit_p[state][obs[n]]+trans_p[state][pState] + uk[n-1][state]
    #             if currProb > currMaxProb:
    #                 currMaxProb = currProb
    #                 prevSelState = pState
    #         max_prob = currMaxProb + emit_p[state][obs[n]]
    #         uk[n][state] = max_prob
    #         ukprev[n][state] = prevSelState
    # z = max(uk[-1],key=uk[-1].get)
    # prev = ukprev[-1][z]
    # backtrack = [z,prev]
    # for n in range(N-2,-1,-1):
    #     z = max(uk[-n],key=uk[-n].get)
    #     prev = ukprev[-n][z]
    #     backtrack.append(prev)
    # backtrack = backtrack[::-1]
    return obs, backtrack
