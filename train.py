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
import pickle
import numpy as np
from scipy.stats import geom
import pickle
from itertools import product
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy import stats



def func(l,E,k):
    kmers = [''.join(group) for group in product('ATCG',repeat=k)]
    val = 1/(4**k)
    summation = 0
    for i in kmers:
        summation+=np.e**((E['+'][i]-E['-'][i])*l)
    return (val *summation)-1
def bisection(f,a,b,N,E,k):
    if f(a,E,k)*f(b,E,k) >= 0:
        return None
    a_n = a
    b_n = b
    for n in range(1,N+1):
        m_n = (a_n + b_n)/2
        f_m_n = f(m_n,E,k)
        if f(a_n,E,k)*f_m_n < 0:
            a_n = a_n
            b_n = m_n
        elif f(b_n,E,k)*f_m_n < 0:
            a_n = m_n
            b_n = b_n
        elif f_m_n == 0:
            print("Found exact solution.")
            return m_n
        else:
            print("Bisection method fails.")
            return None
    return (a_n + b_n)/2
def bitscore(score,length,l,c):
    return ((l*score) - (c*np.log(length))/np.log(2))

# Initialize program arguments, see help= for explanation of each
parser = argparse.ArgumentParser()
parser.add_argument("--query",type=str,help='Path to kmer count file for sequences of interest (e.g. functional regions of a ncRNA)')
parser.add_argument('--null', type=str,help='Path to kmer count file that compose null model (e.g. transcriptome, genome, etc.)')
parser.add_argument('--qT',type=float,help='Probability of query to query transition',default=.999)
parser.add_argument('--nT',type=float,help='Probability of null to null transition',default=.9999)
parser.add_argument('--qPrefix',type=str,help='String, Output file prefix;default=None',default='query')
parser.add_argument('--nPrefix',type=str,help='String, Output file prefix;default=None',default='null')
parser.add_argument('--dir',type=str,help='Output directory',default='./')
parser.add_argument('-k',type=int,help='Value of k',default='2')
parser.add_argument('-a',type=str,help='String, Alphabet to generate k-mers (e.g. ATCG); default=ATCG',default='ATCG')
args = parser.parse_args()

# Check if specified directory exists
# If yes, prompt if replace, else end program
# If no, crash
# Else, loop
if not args.dir.endswith('/'):
    args.dir+='/'
newDir = f'{args.dir}{args.qPrefix}_{args.nPrefix}/'
if not os.path.exists(newDir):
    os.mkdir(newDir)
else:
    flag = True
    while flag:
        usrIN = input(f'Directory {newDir} exists, continue? y/n: ').strip().lower()
        if usrIN == 'y':
            flag = False
        elif usrIN == 'n':
            print('Initiating self-destruct sequence')
            sys.exit()
        else:
            print('Please enter y or n')

a,b = .0000001,1
alphabet = [letter for letter in args.a]
k = args.k
# Load k-mer counts
qCount = pickle.load(open(args.query,'rb'))
nCount = pickle.load(open(args.null,'rb'))

# Loop through specified values of k
# Check if they exist in the counts file,
# and call corefunctions.HMM to generate the HMM matrices
if (k in qCount.keys()) and (k in nCount.keys()):
    qKCount = qCount[k]
    nKCount = nCount[k]
    kDir = newDir+f'{k}/'
    if not os.path.exists(kDir):
        os.mkdir(kDir)
    A,E,states,pi = corefunctions.HMM(qKCount,nKCount,k,args.a,args.qT,args.nT)
    kmers = [''.join(p) for p in itertools.product(alphabet,repeat=k)]

    llr = [E['+'][i]-E['-'][i] for i in kmers]
    llr = np.array(llr)
    l = bisection(func,a,b,1000,E,k)
    adjaLLR = llr*l


    haha = defaultdict(list)
    for i in range(1,100):
        randL = i
        for i in range(10000):
            score = np.zeros(randL)
            scores = np.random.choice(adjaLLR,randL,)
            scores = np.sum(scores)
            if scores>0:
                haha[randL].append(scores)

    means = np.zeros(len(haha))
    C = []
    for j,i in enumerate(haha):
        data = np.array(haha[i])
        mean = np.mean(data)
        means[j] = mean
    logfit = np.polyfit(np.log(np.linspace(1,99,99)), means,1)
    haha = []

    for i in range(1000000):
        randL = geom.rvs(p=.01)
        score = np.zeros(randL)
        scores = np.random.choice(llr,randL,)
        scores = np.sum(scores)
        if scores>0:
            haha.append((scores,randL))
    test = []
    c = logfit[0]
    for i in range(len(haha)):
        test.append(bitscore(haha[i][0],haha[i][1],l,c))
    test = np.array(test)
    scores = test[test>0]
    fit = stats.genpareto.fit(scores)
        # queryMkv = corefunctions.transitionMatrix(qKCount,k,alphabet)
    # nullMkv = corefunctions.transitionMatrix(nKCount,k,alphabet)
    # lgTbl = corefunctions.logLTbl(queryMkv,nullMkv)
else:
    print(f'Missing {k}-mer counts in count file... skipping')

    # np.savetxt(f'{kDir}logtbl.mkv',lgTbl)
pickle.dump({'A':A,'E':E,'pi':pi,'states':states,'gp':fit,'l':l,'c':logfit[0]},open(f'{kDir}hmm.mkv','wb'))
