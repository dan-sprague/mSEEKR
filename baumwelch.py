from seekr.fasta_reader import Reader
import numpy as np
import corefunctions
import argparse
import pickle
from tqdm import tqdm as tqdm
from scipy.special import logsumexp
from collections import defaultdict


parser = argparse.ArgumentParser()
parser.add_argument("-k",type=int)
parser.add_argument('--db',type=str,help='Path to fasta file with sequences to calculate similarity score')
parser.add_argument('--prior',type=str)
args = parser.parse_args()


fa = Reader(args.db)
seqs = fa.get_seqs()[0]
model = args.prior

k = args.k
O,oIdx,nBP = corefunctions.kmersWithAmbigIndex(seqs,k)
hmm = pickle.load(open(args.prior,'rb'))
A,E,pi,states = hmm['A'],hmm['E'],hmm['pi'],hmm['states']

data = defaultdict(list)
for i in tqdm(range(100)):
    a = corefunctions.fwd(O,A,pi,states,E)
    b = corefunctions.bkw(O,A,pi,states,E)
    A = corefunctions.update(a,b,O,states,A,E)
    for i in states:
        marg = logsumexp(list(A[i].values()))
        A[i]['+']-=marg
        A[i]['-']-=marg
    data['alpha'].append(A['+']['+'])
    data['beta'].append(A['-']['-'])

pickle.dump(data,open(f'./{args.prior}_{k}_bwMLE.p','wb'))