import kmers
import argparse
import pickle
import os
from multiprocessing import pool
from itertools import product
from seekr.fasta_reader import Reader

parser = argparse.ArgumentParser()
parser.add_argument('--fasta', type=str,help='Path to fasta file')
parser.add_argument('--name',type=str,help='name for count file')
parser.add_argument('-k',type=str,help='Comma delimited string of possible k-mer values',default='2,3,4')
parser.add_argument('-a',type=str,help='String, Alphabet to generate k-mers (e.g. ATCG); default=ATCG',default='ATCG')
parser.add_argument('-n',type=int,help='Number of CPU cores. Each job corresponds to a value of k, and the program scales well with multiprocessing',default=1)

args = parser.parse_args()
kVals = [int(i) for i in args.k.split(',')]
a = args.a.upper()

F = Reader(args.fasta)
fS = F.get_seqs()

fString = '$'.join(fS)

with pool.Pool(args.n) as multiN:
    jobs = multiN.starmap(kmers.main,product(*[[fString],kVals,[a]]))
    dataDict = dict(jobs)
kDir = './counts/'
pickle.dump(dataDict,open(f'{kDir}{args.name}.skr','wb'))
