import corefunctions
import manaStats
import argparse
import itertools
import numpy as np
from seekr.fasta_reader import Reader
from scipy.stats import norm
from collections import defaultdict
from tqdm import tqdm
parser = argparse.ArgumentParser()
parser.add_argument("--query",type=str,help='Target sequences')
parser.add_argument('--null', type=str,help='Sequences that compose null model')
parser.add_argument('--db',type=str,help='Fasta file with sequences to calculate similarity score')
parser.add_argument('-p', type=int,help='Solve for score S that integrates KDE to this p-val. If S < 0 (S is more rare than anticipated), set S = 0',default=.01)
parser.add_argument('-o',type=int,help='Order of markov model',default=3)
parser.add_argument('-w', type=int, help='Window for tile size', default=200)
parser.add_argument('-s', type=int, help='How many bp to slide tiles', default=20)
parser.add_argument('-a',type=str,help='Alphabet to generate k-mers',default='ATCG')
parser.add_argument('--gtf',type=str,help='GTF file of loci of interest',default=None)
parser.add_argument('--narrowPeak',type=str,help='narrowPeak file of protein binding peaks',default=None)
parser.add_argument('--prefix',type=str,help='Prefix to append to filenames',default='out')
parser.add_argument('-n',type=int,help='Number of cores to use, default = 1',default=1)
args = parser.parse_args()

alphabet = [letter for letter in args.a]
kmers = [''.join(p) for p in itertools.product(alphabet,repeat=args.o)]
queryMkv = corefunctions.trainModel(args.query,args.o,kmers,alphabet)
nullMkv = corefunctions.trainModel(args.null,args.o,kmers,alphabet)

lgTbl = corefunctions.logLTbl(queryMkv,nullMkv)
np.savetxt(f'{args.prefix}_{args.o}order_LRTbl.mkv',lgTbl)

probMap = {'A':.3,'T':.3,'C':.2,'G':.2}
probs = [probMap[letter] for letter in args.a]
randSeqs = [manaStats.dnaGen(args.w,alphabet,probs) for i in range(1000)]
randSeqsScore = np.array([corefunctions.classify(seq,args.o,lgTbl,alphabet) for seq in randSeqs])
kde = manaStats.KDE(randSeqsScore.reshape(-1,1))

S = manaStats.kdeCDF(kde.best_estimator_,1000,-100,100,args.p)
if S < 0:
    S = 0

target = Reader(args.db)
targetSeqs,targetHeaders = target.get_seqs(),target.get_headers()

targetMap = defaultdict(list)
for tHead,tSeq in zip(targetHeaders,targetSeqs):
    tileScores = np.array([corefunctions.classify(tSeq[i:i+args.w],args.o,lgTbl,alphabet) for i in range(0,len(tSeq),args.s)])
    randSums = np.zeros(30)
    for i in range(30):
        samp = np.array(kde.best_estimator_.sample(len(tileScores)))
        randSums[i] = np.sum(samp[samp>0])
    normSums = norm(np.mean(randSums),np.std(randSums))
    P = 1-normSums.cdf(np.sum(tileScores[tileScores>0]))

    targetMap[tHead].append([np.sum(tileScores[tileScores>0]),P,np.sum(tileScores>0),manaStats.tileE(tileScores,args.p,np.sum(tileScores>0))])

    for i,tileScore in tqdm(enumerate(tileScores)):
        if tileScore > S:
            targetMap[tHead].append(f'{i}\t{i*args.s}:{(i*args.s)+args.w}\t{tSeq[i*args.s:(i*args.s)+args.w]}\t{tileScore}\t{manaStats.integrate(kde.best_estimator_,1000,-100,tileScore)[0]}\n')

with open('./align.txt','w') as outfile:
    for h,data in targetMap.items():
        outfile.write(f'$ {tHead}\t{data[0]}\n')
        outfile.write(f'Tile\tbp Range\tSequence\tLog-Likelihood\tp-val\n')
        for string in data[1:]:
            outfile.write(string)
