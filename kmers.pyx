from seekr.fasta_reader import Reader
from itertools import product


cpdef main(str fString,int k,str a):
    cdef int K = k
    cdef list kmers = [''.join(p) for p in product(a,repeat=K)]
    cdef int pos = 4**K
    cdef int j
    cdef dict kmerDict = {}
    for j in range(pos):
        kmerDict[kmers[j]] = 1
    cdef int N = len(fString)
    r = countKmers(fString,N,K,kmerDict)
    return k,r

cdef countKmers(str fString,int N,int K,dict kmerDict):
    cdef int i
    for i in range(N-K+1):
        curr = fString[i:i+K]
        if curr in kmerDict:
            kmerDict[curr]+=1
    return kmerDict
