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
import matplotlib.pyplot as plt
import seaborn as sns
parser = argparse.ArgumentParser()
parser.add_argument('--mkv',type=str,help='Directory containing markov (.mkv) models',default='./markovModels/')
parser.add_argument('--db',type=str,help='Path to fasta file with sequences to calculate similarity score')
parser.add_argument('--win',type=str,help='Comma delimited string (no quotes) of window sizes; default=100,200',default='100,200')
parser.add_argument('--slide',type=str,help='Comma delimited string (no quotes) of how much to slide for each window size specified. If a single number is specified, that slide will be applied to all windows, otherwise, the length of this list must equal the window argument length',default='20')
parser.add_argument('--dir',type=str,help='Path to save directory; default=./',default='./pdfs/')
parser.add_argument('-a',type=str,help='Alphabet, string (no quotes); default=ATCG',default='ATCG')
parser.add_argument('-e',type=str,help='Experiment name (no quotes)')
args = parser.parse_args()

w = [int(i) for i in args.win.split(',')]
slide = args.slide.split(',')
if len(slide)==1:
    slide = slide*len(w)
slide = [int(i) for i in slide]
windowSlide = list(zip(w,slide))
alphabet = [letter for letter in args.a]

models = [f for f in glob.iglob(f'{args.dir}*.mkv')]
models = sorted(models)
modelDict = defaultdict(list)
for model in models:
    prefix = os.path.basename(model).split('_')[:2]
    query = prefix[0]
    null = prefix[1]
    modelDict[query+null].append(model)
dataBase = Reader(args.db)
dataBase = zip(dataBase.get_headers(),dataBase.get_seqs())
sns.set_context('talk',font_scale=.5)
for tH,tSeq in dataBase:
    for modelGroup,models in modelDict.items():
        fig,axes = plt.subplots(nrows=len(w),ncols=len(models),figsize=(7,12))
        for col,model in enumerate(models):
            lgTbl = np.loadtxt(model)
            k = int(log(lgTbl.size,len(args.a)))
            prefix = os.path.basename(model).split('_')[:2]
            for row,(currWin,currSlide) in enumerate(windowSlide):
                tiles = [tSeq[i:i+currWin] for i in range(0,len(tSeq)-currWin+1,currSlide)]
                tileScores = np.array([corefunctions.score(tile,k,lgTbl,alphabet) for tile in tiles])
                axes[row,col].plot(tileScores)
                axes[row,col].axhline(y=0,linestyle='--')

        fig.tight_layout()
        fig.savefig(f'{args.e}_{tH[:5]}_{modelGroup}.pdf')
        plt.close()
    # plotIdx =0
    # for row in ax:
    #     for col in row:
    #         col.plot(modelData[plotIdx])
    #         col.axhline(y=0,linestyle='--')
    #         plotIdx+=1
    # plt.tight_layout()
    # plt.savefig(f'{args.dir}{prefix}_{args.e}.pdf',bbox_inches='tight')
    # plt.close()
