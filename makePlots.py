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

parser = argparse.ArgumentParser()

parser.add_argument('--db',type=str,help='Path to fasta file with sequences to calculate similarity score')
parser.add_argument('--dir',type=str,help='Path to save directory; default=./',default='./')
parser.add_argument('-a',type=str,help='Alphabet; default=ATCG',default='ATCG')
parser.add_argument('-e',type=str,help='Experiment name')
args = parser.parse_args()

w = [(25,20),(50,20),(100,20),(200,20),(500,20)]
alphabet = [letter for letter in args.a]

models = [f for f in glob.iglob('./*mkv')]
for model in models:
    for win,sl in w:
        lgTbl = np.loadtxt(model)
        # Explicitly determine k from the size of the log matrix and the size of the alphabet used to generate it
        k = int(log(lgTbl.size,len(args.a)))
        prefix = os.path.basename(model).split('_')[:2]
        dataBase = Reader(args.db)
        dataBase = zip(dataBase.get_headers(),dataBase.get_seqs())
        for tH,tSeq in dataBase:
            tiles = [tSeq[i:i+win] for i in range(0,len(tSeq)-win+1,sl)]
            tileScores = np.array([corefunctions.score(tile,k,lgTbl,alphabet) for tile in tiles])
            corefunctions.plotTiles(tileScores,f'/mnt/c/Users/sprag/Documents/{args.e}_{prefix}_{k}_{win}w_{sl}sl_tilePlot.pdf',0)
            plt.close('all')
