import corefunctions
import numpy as np
from seekr.fasta_reader import Reader
import argparse
import os
parser = argparse.ArgumentParser()
parser.add_argument("--model",type=str,help='Path to directory containing .mkv files or path to a single .mkv file;default=./markovModels/',default='./markovModels/')
parser.add_argument('--null',type=str,help='Path to fasta file with sequences to calculate similarity score')
parser.add_argument('--nullPrefix',type=str,default='null')
parser.add_argument('--modelPrefix',type=str,default='model')
parser.add_argument('--queryPrefix',type=str,default='query')
parser.add_argument('-q',type=str,help='query')
parser.add_argument('-k',type=int,help='k')
args = parser.parse_args()

query = Reader(args.q)
queryHeader,querySeq = query.get_headers(),query.get_seqs()

S = [args.model,args.null]

alphabet = ['A','T','C','G']

seq = []

A,E,states,pi = corefunctions.HMM(S,args.k,alphabet)

with open(f'./{args.queryPrefix}_{args.modelPrefix}_{args.nullPrefix}_{args.k}.txt','w') as outfile:
    for qH,obs in zip(queryHeader,querySeq):
        O = [obs[i:i+args.k] for i in range(0,len(obs)-args.k+1)]
        bTrack = corefunctions.viterbi(O,A,E,states,pi)
        currSeq = ''
        for o,s in zip(O,bTrack):
            if s == 'Model':
                currSeq+='1'
            elif s=='Null':
                currSeq+='0'
        outfile.write(f'$ {qH}\n')
        outfile.write(currSeq)


# import matplotlib.pyplot as plt
# plt.figure()
# for i,s in enumerate(seq):
#     if s == '+':
#         plt.scatter(x=i,y=0,marker='s',color='r')
#     if s=='-':
#         plt.scatter(x=i,y=0,marker='s',color='blue')
# plt.show()
