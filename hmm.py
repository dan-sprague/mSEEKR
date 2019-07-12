import corefunctions
import numpy as np
from seekr.fasta_reader import Reader
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--model",type=str,help='Path to directory containing .mkv files or path to a single .mkv file;default=./markovModels/',default='./markovModels/')
parser.add_argument('--null',type=str,help='Path to fasta file with sequences to calculate similarity score')
parser.add_argument('-q',type=str,help='query')
parser.add_argument('-k',type=int,help='k')
args = parser.parse_args()

query = Reader(args.q)
querySeq = ''.join(query.get_seqs())

S = [args.model,args.null]

alphabet = ['A','T','C','G']


obs,states,start_p,trans_p,emit_p = corefunctions.HMM(S,args.k,alphabet,querySeq)
bTrack = corefunctions.viterbi(obs,states,start_p,trans_p,emit_p)
seq = ''

for o,s in zip(obs,bTrack):
    if s == 'Model':
        seq+='+'
    elif s=='Null':
        seq+='-'
print(seq)
import matplotlib.pyplot as plt
plt.figure()
for i,s in enumerate(seq):
    if s == '+':
        plt.scatter(x=i,y=0,marker='s',color='r')
    if s=='-':
        plt.scatter(x=i,y=0,marker='s',color='blue')
plt.show()
