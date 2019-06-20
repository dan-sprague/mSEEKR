from scipy.stats import poisson
from sklearn.neighbors import KernelDensity
#from sklearn.model_selection import GridSearchCV
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
from scipy.integrate import simps
from scipy.stats import gaussian_kde
def dnaGen(length,alphabet,probs):
    DNA = choice(alphabet,length,p=probs)
    return ''.join(DNA)

# def slvLimit(K,n,a,b,p):
#     h = (b-a)/n
#     A = .5 * h * (f(a,K))
#     count =1
#     while A <= (1-p):
#         A += h * f(a+(count*h),K)
#         count+=1
#     return a + (count*h)
# def integrate(b,x,y):
#     x = x[x<=b]
#     y = y[:len(x)]
#     AUC = simps(y,x)
#     return 1.0-AUC
def tileE(x,P,O):
    E = len(x)*P
    rv = poisson(E)
    return 1-rv.cdf(O)
# def findS(x,y,p):
#     idx = np.abs(y-p).argmin()
#     return x[idx]
