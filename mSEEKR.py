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

def calculateSimilarity(data):
    tHead,tSeq = data
    tiles = [tSeq[i:i+args.w] for i in range(0,len(tSeq)-args.w+1,args.s)]
    tileScores = np.array([corefunctions.score(tile,k,lgTbl,alphabet) for tile in tiles])
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
    bpHits = []
    for i in argSortScores:
        bpHits.append(list(range((i*args.s),(i*args.s)+args.w)))
        tileScore = tileScores[i]
        # integratedP = kde.integrate_box(tileScore,upperLimit*10)
        #integratedP = 1-testNorm.cdf(tileScore)
        str1 = f'{i}\t{i*args.s}:{(i*args.s)+args.w}\t'
        str2 = f'{tSeq[i*args.s:(i*args.s)+args.w]}\t{tileScore}\n'
        # str3 = f'{integratedP}\n'
        strData = str1+str2
        strDataList.append(strData)
    bpHits = np.sort(np.unique(np.array(bpHits).flatten()))
    return tHead,[summaryStats,strDataList,[''.join([tSeq[i] for i in bpHits])]]


parser = argparse.ArgumentParser()
parser.add_argument("--model",type=str,help='Path to directory containing .mkv files or path to a single .mkv file;default=./markovModels/',default='./markovModels/')
parser.add_argument('--db',type=str,help='Path to fasta file with sequences to calculate similarity score')
parser.add_argument('--nRAND',type=int,help='Int >0, Number of random sequences to generate. default=10^5',default=10**5)
parser.add_argument('--prefix',type=str,help='String, Output file prefix;default=None')
parser.add_argument('--retrain',help='Save unique sequence hits to fasta file for markov training',action='store_true')
parser.add_argument('--bkg',type=str,help='Path to fasta file from which to calculate background nucleotide frequencies, if not passed default is uniform',default=None)
parser.add_argument('-p', type=float,help='Float, Desired p-val of log-likelihood ratio score "S" to consider significant. If P(S > 0) is very small and less than this argument, set S = 0 and p = P(S>0); default=.01',default=.01)
parser.add_argument('-w', type=int, help='Integer >= k, length of tile sizes; default=200', default=200)
parser.add_argument('-s', type=int, help='Integer >=1, how many bp to slide tiles. Increasing this parameter decreases compute time significantly; default=20', default=20)
parser.add_argument('-a',type=str,help='String, Alphabet to generate k-mers (e.g. ATCG); default=ATCG',default='ATCG')
parser.add_argument('-n',type=int,help='Integer 1 <= n <= max(cores), Number of processor cores to use; default = 1. This scales with the number of sequence comparisons in --db',default=1)


args = parser.parse_args()
alphabet = [letter for letter in args.a]

if os.path.isdir(args.model):
    models = [f for f in glob.iglob(f'{args.model}*mkv')]
else:
    models = [args.model]

if args.bkg:
    bkgFa = Reader(args.bkg)
    bkgSeqs = bkgFa.get_seqs()
    probMap = corefunctions.nucContent(bkgSeqs,args.a)
    print(f'Background Frequencies: {probMap}')
elif not args.bkg:
    probMap = {'A':.25,'T':.25,'C':.25,'G':.25}


for model in models:
    modelName = os.path.basename(model)
    modelName = '_'.join(modelName.split('_')[:2])
    lgTbl = np.loadtxt(model)
    # Explicitly determine k from the size of the log matrix and the size of the alphabet used to generate it
    k = int(log(lgTbl.size,len(args.a)))
    kmers = [''.join(p) for p in itertools.product(alphabet,repeat=k)]
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
    #minP = 1-kde.integrate_box_1d(lowerLimit,upperLimit)
    #minPStr = f'< 1E{np.ceil(np.log10(abs(minP)))}'
    S = x[i]
    print(f'Score Threshold: {S}\nEst. p-val: {1-kde.integrate_box(lowerLimit,S)}')
    # If P(S > 0) < args.p, set S = 0
    if S < 0:
        S = 0
        args.p = kde.integrate_box(S,upperLimit*10)
        print(f'S < 0, setting S = 0\np-val: {args.p}')
    print('\nDone')
    target = Reader(args.db)
    targetSeqs,targetHeaders = target.get_seqs(),target.get_headers()

    targetMap = defaultdict(list)
    print('\nScanning database sequences')
    with pool.Pool(args.n) as multiN:
        jobs = multiN.starmap(calculateSimilarity,product(*[list(zip(targetHeaders,targetSeqs))]))
        dataDict = dict(jobs)
    print('\nDone')
    with open(f'./{args.prefix}_{modelName}_{k}_{args.w}w_{args.s}sl_HSS.txt','w') as outfile:
        for h,data in dataDict.items():
            outfile.write(f'$ {h}\t{data[0]}\n')
            outfile.write(f'Tile\tbp Range\tSequence\tLog-Likelihood\tp-val\n')
            for string in data[1]:
                outfile.write(string)
    if args.retrain:
        with open(f'./{args.prefix}_{modelName}_{k}_{args.w}w_{args.s}sl.fa','w') as outfile:
            for h, data in dataDict.items():
                outfile.write(f'>{h}\n')
                for seq in data[2]:
                    outfile.write(f'{seq}\n')
