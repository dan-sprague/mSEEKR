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

def mainCompute(data):
    tHead,tSeq = data
    tiles = [tSeq[i:i+args.w] for i in range(0,len(tSeq)-args.w+1,args.s)]
    tileScores = np.array([corefunctions.score(tile,k,lgTbl,alphabet) for tile in tiles])
    corefunctions.plotTiles(tileScores,f'/mnt/c/Users/sprag/Documents/{args.prefix}_{k}_{args.w}w_{args.s}sl_tilePlot.pdf',S)
    randSums = np.zeros(args.nRAND)
    for i in range(args.nRAND):
        samp = np.array(kde.resample(len(tileScores)))
        randSums[i] = np.sum(samp[samp>S])
    hitSum = np.sum(tileScores[tileScores>S])
    if any(randSums>0):
        normSums = norm(np.mean(randSums),np.std(randSums))
        sumP = 1-normSums.cdf(hitSum)
    # No sums greater than S were observed
    else:
        sumP=0
    summaryStats = [hitSum,sumP,np.sum(tileScores>S),coreStats.tileE(tileScores,args.p,np.sum(tileScores>S))]
    argSortScores = np.argsort(tileScores)[::-1]
    idxHit = np.nonzero(tileScores>S)
    argSortScores = argSortScores[np.isin(argSortScores,idxHit)]
    strDataList = []
    for i in argSortScores:
        tileScore = tileScores[i]
        integratedP = 1-kde.integrate_box(lowerLimit,tileScore)
        if integratedP <= 0:
            integratedP = minPStr
        #integratedP = 1-testNorm.cdf(tileScore)
        str1 = f'{i}\t{i*args.s}:{(i*args.s)+args.w}\t'
        str2 = f'{tSeq[i*args.s:(i*args.s)+args.w]}\t{tileScore}\t'
        str3 = f'{integratedP}\n'
        strData = str1+str2+str3
        strDataList.append(strData)
    return tHead,[summaryStats,strDataList]

parser = argparse.ArgumentParser()
parser.add_argument("--model",type=str,help='Path to fasta file or containing sequences to build markov model (e.g. functional regions of a ncRNA)')
parser.add_argument('--db',type=str,help='Path to fasta file with sequences to calculate similarity score')
parser.add_argument('--nRAND',type=int,help='Int >0, Number of random sequences to generate. If using empircal p-vals, minimum p-val is 1/nRAND; default=10^5',default=10**5)
parser.add_argument('--prefix',type=str,help='String, Output file prefix;default=None')
parser.add_argument('-p', type=int,help='Float, Desired p-val of log-likelihood ratio score "S" to consider significant. If P(S > 0) is very small and less than this argument, set S = 0 and p = P(S>0); default=.01',default=.01)
parser.add_argument('-w', type=int, help='Integer >= k, length of tile sizes; default=200', default=200)
parser.add_argument('-s', type=int, help='Integer >=1, how many bp to slide tiles. Increasing this parameter decreases compute time significantly; default=20', default=20)
parser.add_argument('-a',type=str,help='String, Alphabet to generate k-mers (e.g. ATCG); default=ATCG',default='ATCG')
parser.add_argument('-n',type=int,help='Integer 1 <= n <= max(cores), Number of processor cores to use; default = 1',default=1)


args = parser.parse_args()

lgTbl = np.loadtxt(args.model)
# Explicitly determine k from the size of the log matrix and the size of the alphabet used to generate it
k = int(log(lgTbl.size,len(args.a)))
alphabet = [letter for letter in args.a]
kmers = [''.join(p) for p in itertools.product(alphabet,repeat=k)]
probMap = {'A':.3,'T':.3,'C':.2,'G':.2}
probs = [probMap[letter] for letter in args.a]
print('\nGenerating model of score distribution')
randSeqs = [coreStats.dnaGen(args.w,alphabet,probs) for i in range(args.nRAND)]
randSeqsScore = np.array([corefunctions.score(seq,k,lgTbl,alphabet) for seq in randSeqs])
kde = gaussian_kde(randSeqsScore)
lowerLimit= np.min(lgTbl) * args.w
upperLimit = np.max(lgTbl) * args.w
x = np.linspace(lowerLimit,upperLimit,10000)
F = 0
i=1
while F < (1 - args.p):
    F=kde.integrate_box(lowerLimit,x[i])
    i+=1
minP = 1-kde.integrate_box_1d(lowerLimit,upperLimit)
minPStr = f'< 1E{np.ceil(np.log10(abs(minP)))}'
S = x[i]
print(f'Score Threshold: {S}\nEst. p-val: {1-kde.integrate_box(lowerLimit,S)}')
# If P(S > 0) < args.p, set S = 0
if S < 0:
    S = 0
    args.p = 1-kde.integrate_box(lowerLimit,S)
    if args.p <= 0:
        args.p = 10**(np.ceil(np.log10(abs(minP))))
    print(f'S < 0, setting S = 0\np-val: {args.p}')
print('\nDone')
target = Reader(args.db)
targetSeqs,targetHeaders = target.get_seqs(),target.get_headers()

targetMap = defaultdict(list)
print('\nScanning database sequences')
with pool.Pool(args.n) as multiN:
    jobs = multiN.starmap(mainCompute,product(*[list(zip(targetHeaders,targetSeqs))]))
    dataDict = dict(jobs)
print('\nDone')
with open(f'./{args.prefix}_{k}_{args.w}w_{args.s}sl_HSS.txt','w') as outfile:
    for h,data in dataDict.items():
        outfile.write(f'$ {h}\t{data[0]}\n')
        outfile.write(f'Tile\tbp Range\tSequence\tLog-Likelihood\tp-val\n')
        for string in data[1]:
            outfile.write(string)
