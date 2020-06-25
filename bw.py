
from seekr.fasta_reader import Reader
import numpy as np
import corefunctions
import argparse
import pickle
from tqdm import tqdm as tqdm
from scipy.special import logsumexp
from collections import defaultdict


```
Author: Daniel Sprague
Calabrese Lab, Department of Pharmacology UNC Chapel Hill
June 24, 2020

Purpose: Run the Baum-Welch MLE algorithm to identify putative transition parameters for mSEEKR given a set of training sequences.
    
Input: 
    
-k : K-mer value to use
--db : Training sequences
--prior : Output of train.py in the mSEEKR pipeline, essentially, initial guesses for the transition parameters, 
    as well as trained values for the k-mer frequencies for query and null states
--its : How many iterations of the BW algorithm to run (default=100)
    
Output:
    
    A new .mkv python pickle file (binary) that contains the original k-mer frequencies and updated transition matrix 

```


# Initiate command line arguments


parser = argparse.ArgumentParser()
parser.add_argument("-k",type=int)
parser.add_argument('--db',type=str,help='Path to fasta file containing training sequences')
parser.add_argument('--prior',type=str)
parser.add_argument('--its',type=int,help='Iterations to do, default=100',default=100)
args = parser.parse_args()


fa = Reader(args.db)
seqs = fa.get_seqs()[0]
model = args.prior

k = args.k

# Identify the location of any ambiguous nucleotides (N)
O,oIdx,nBP = corefunctions.kmersWithAmbigIndex(seqs,k)
# Load in train.py output 
hmm = pickle.load(open(args.prior,'rb'))

```
A - transition matrix (dictionary)
E - emission matrix (dictionary)
pi - initial state probabilities (always 50/50)
states - list of states (query,null)
```
A,E,pi,states = hmm['A'],hmm['E'],hmm['pi'],hmm['states']

data = defaultdict(list)

# Perform --its number of baum-welch iterations
# BW converges to a LOCAL MINIMUM -- best to perform a grid search (run this script on a variety of different --priors)
for i in tqdm(range(args.its)):
    a = corefunctions.fwd(O,A,pi,states,E) # forward probabilities
    b = corefunctions.bkw(O,A,pi,states,E) # backwards probabilities
    A = corefunctions.update(a,b,O,states,A,E) # update parameter values
    
    # The above step yields probabilities that are not yet normalized (do not sum to 1)
    # That is, for prior state "i" (can be either query or null) transitioning to either query (+) or null (-): p(-|i) + p(+|j) != 1 
    # Here, we calculate the marginal probability for transitioning from state i to either state i or j (normalizing constant)
    # Then, divide each probability p(-|i) and p(+|i) by the marginal probability 
    # Now, p(-|i) + p(+|j) = 1
    for i in states:
        # marginal probability (normalizing constant)
        marg = logsumexp(list(A[i].values())) # cannot simply sum in log-space, logsumexp to sum log-space numbers without over/under-flow (see wiki page)
        A[i]['+']-=marg # log-space, so this represents the division by the normalizing constant
        A[i]['-']-=marg # as above
    data['alpha'].append(A['+']['+'])
    data['beta'].append(A['-']['-'])
  
# Data dictionary tracks the iterations of BW algorithm for manual inspection if necessary 

pickle.dump(data,open(f'./{args.prior}_{k}_bwMLE.p','wb'))
