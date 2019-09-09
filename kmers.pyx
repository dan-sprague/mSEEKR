from seekr.fasta_reader import Reader
from itertools import product

cpdef main(str fString,int k,str a):
    cdef unsigned int K = k
    cdef list kmers = [''.join(p) for p in product(a,repeat=K)]
    cdef unsigned int pos = 4**K
    cdef unsigned int j
    cdef dict kmerDict = {}
    for j in range(pos):
        kmerDict[kmers[j]] = 1
    cdef unsigned long long N = len(fString)
    r = countKmers(fString,N,K,kmerDict)
    return k,r

cdef countKmers(str fString,unsigned long long N,unsigned int K,dict kmerDict):
    cdef unsigned long long i
    for i in range(N-K+1):
        curr = fString[i:i+K]
        if curr in kmerDict:
            kmerDict[curr]+=1
    return kmerDict
