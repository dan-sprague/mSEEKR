import corefunctions
import manaStats
import argparse
import itertools
import numpy as np
parser = argparse.ArgumentParser()
parser.add_argument("--query",type=str,help='Target sequences')
parser.add_argument('--null', type=str,help='Sequences that compose null model')
parser.add_argument('--db',type=str,help='Fasta file with sequences to calculate similarity score')
parser.add_argument('-p', type=int,help='p-val threshold to consider score "S" significant. Note: If p(S>0) < p, threshold will be set to a score of 0',default=.01)
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
