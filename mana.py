import corefunctions
import manaStats
import argparse
import itertools
import numpy as np
from seekr.fasta_reader import Reader
from scipy.stats import norm
from collections import defaultdict
from tqdm import tqdm
from multiprocessing import pool
parser = argparse.ArgumentParser()
parser.add_argument("--query",type=str,help='Path to fasta file containing sequences to build markov model')
parser.add_argument('--null', type=str,help='Path to fasta file containing sequences that compose null model')
parser.add_argument('--db',type=str,help='Path to fasta file with sequences to calculate similarity score')
parser.add_argument('-p', type=int,help='Desired p-val of score to consider significant. If P(S > 0) is very small and less than this argument, set S = 0 and P = P(S>0)',default=.01)
parser.add_argument('-o',type=int,help='Order of markov model',default=3)
parser.add_argument('-w', type=int, help='Window for tile size', default=200)
parser.add_argument('-s', type=int, help='How many bp to slide tiles', default=20)
parser.add_argument('-a',type=str,help='Alphabet to generate k-mers',default='ATCG')
parser.add_argument('-n',type=int,help='Number of cores to use, default = 1',default=1)
parser.add_argument('--prefix',type=str,help='Output prefix')
args = parser.parse_args()

alphabet = [letter for letter in args.a]
print('Counting k-mers...')
kmers = [''.join(p) for p in itertools.product(alphabet,repeat=args.o)]
queryMkv = corefunctions.trainModel(args.query,args.o,kmers,alphabet)
nullMkv = corefunctions.trainModel(args.null,args.o,kmers,alphabet)
lgTbl = corefunctions.logLTbl(queryMkv,nullMkv)
probMap = {'A':.3,'T':.3,'C':.2,'G':.2}
probs = [probMap[letter] for letter in args.a]
print('Done')
print('\nGenerating model of score distribution')
randSeqs = [manaStats.dnaGen(args.w,alphabet,probs) for i in range(5000)]
randSeqsScore = np.array([corefunctions.classify(seq,args.o,lgTbl,alphabet) for seq in randSeqs])
kde = manaStats.KDE(randSeqsScore.reshape(-1,1))

# Calculate lower limit of integration and integrate KDE
lowerLimit= np.min(lgTbl) * args.w
S = manaStats.kdeCDF(kde,5000,lowerLimit,100,args.p)
print(f'Score Threshold: {S}')
# If P(S > 0) < args.p, set S = 0
if S < 0:
    S = 0
    args.p = manaStats.integrate(kde,5000,lowerLimit,S)[0]
print('\nDone')
target = Reader(args.db)
targetSeqs,targetHeaders = target.get_seqs(),target.get_headers()

targetMap = defaultdict(list)
print('\nScanning database sequences')
for tHead,tSeq in zip(targetHeaders,targetSeqs):
    tileScores = np.array([corefunctions.classify(tSeq[i:i+args.w],args.o,lgTbl,alphabet) for i in range(0,len(tSeq),args.s)])
    randSums = np.zeros(1000)
    for i in range(1000):
        samp = np.array(kde.sample(len(tileScores)))
        randSums[i] = np.sum(samp[samp>S])
    hitSum = np.sum(tileScores[tileScores>S])
    if any(randSums>0):
        normSums = norm(np.mean(randSums),np.std(randSums))
        sumP = 1-normSums.cdf(hitSum)
    # No sums greater than S were observed
    else:
        sumP=0
    targetMap[tHead].append([hitSum,sumP,np.sum(tileScores>S),manaStats.tileE(tileScores,args.p,np.sum(tileScores>S))])
    argSortScores = np.argsort(tileScores)[::-1]
    idxHit = np.nonzero(tileScores>S)
    argSortScores = argSortScores[np.isin(argSortScores,idxHit)]
    for i in argSortScores:
        tileScore = tileScores[i]
        integratedP = manaStats.integrate(kde,5000,lowerLimit,tileScore)[0]
        str1 = f'{i}\t{i*args.s}:{(i*args.s)+args.w}\t'
        str2 = f'{tSeq[i*args.s:(i*args.s)+args.w]}\t{tileScore}\t'
        str3 = f'{integratedP}\n'
        strData = str1+str2+str3
        targetMap[tHead].append(strData)
print('\nDone')
with open(f'./{args.prefix}_{args.o}o_{args.w}w_{args.s}sl_HSS.txt','w') as outfile:
    for h,data in targetMap.items():
        outfile.write(f'$ {tHead}\t{data[0]}\n')
        outfile.write(f'Tile\tbp Range\tSequence\tLog-Likelihood\tp-val\n')
        for string in data[1:]:
            outfile.write(string)
