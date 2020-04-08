import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from seekr.kmer_counts import BasicCounter
from scipy.stats import entropy

def rectCorr(S,Q):

    S = S-np.mean(S,axis=1,keepdims=True)
    Q = Q-np.mean(Q,axis=1,keepdims=True)
    covariance = S.dot(Q.T)
    sN = np.linalg.norm(S,axis=1)
    qN = np.linalg.norm(Q,axis=1)
    norm = np.outer(sN,qN)

    CORR = covariance/norm
    return CORR.ravel()

d = {}
for k in [2,3,4,5]:
    ref = BasicCounter(infasta='../../mSEEKR/fastaFiles/gencode.vM17.lncRNA_transcripts.fa',k=k,std=False,mean=False,log2=False,silent=True)
    ref.get_counts()
    R = np.log2(ref.counts+1)
    for query in ['mA','mB','hD','mE']:
        queryC = BasicCounter(infasta=f'./fastaFiles/{query}.fa',k=k,std=False,mean=False,log2=False,silent=True)
        queryC.get_counts()
        Q = np.log2(queryC.counts+1)
        Q = (Q-np.mean(R,axis=0))/np.std(R,axis=0)
        rC = rectCorr(R,Q)
        for mK in [2,3,4,5]:
            for p in [5,75,9,99,999,9999]:
                for q in [5,75,9,99,999,9999]:
                    try:
                        df = pd.read_csv(f'./mouse_{p}_{q}_{query}_{p}_{q}_mT_{mK}.txt.skr_{k}.txt',sep='\t',header=None,index_col=0)
                    except:
                        continue
                    print(df.head())
                    df.columns=['start','stop','length','score','Pearson']
                    qH = np.histogram(df['Pearson'].values,bins=20,range=(-1,1))[0]
                    rH = np.histogram(rC,bins=20,range=(-1,1))[0]
                    KL = entropy(qH,rH)
                    print(KL)
                    1/0
