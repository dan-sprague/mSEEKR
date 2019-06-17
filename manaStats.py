from scipy.stats import poisson
from sklearn.neighbors import KernelDensity
from sklearn.model_selection import GridSearchCV
import numpy as np
import itertools
from itertools import product
from seekr.kmer_counts import BasicCounter
from seekr.fasta_reader import Reader
import matplotlib.pyplot as plt
from itertools import groupby
from operator import itemgetter
from collections import defaultdict
from numpy.random import choice
from scipy.stats import norm
from scipy.stats import poisson

def dnaGen(length,alphabet,probs):
    DNA = choice(alphabet,length,p=probs)
    return ''.join(DNA)
def KDE(arr):
    grid = GridSearchCV(KernelDensity(),
                    {'bandwidth': np.linspace(0.1, 3, 30)},
                    cv=10)
    return grid.fit(arr)
def kdeCDF(K,n,a,b,p):
    h = (b-a)/n
    A = .5 * h * (f(a,K)+f(b,K))
    count =0
    while A <= (1-p):
        A += h * f(a+(count*h),K)
        count+=1
    return a + (count*h)
def integrate(K,n,a,b):
    h = (b-a)/n
    A = .5 * h * (f(a,K)+f(b,K))
    for i in range(1,n):
        A += h * f(a+(i*h),K)
    return 1-A
def f(x,kde):
    return np.exp(kde.score_samples(np.array(x).reshape(-1,1)))

def tileE(x,P,O):
    E = len(x)*P
    rv = poisson(E)
    return 1-rv.cdf(O)
