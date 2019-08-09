import itertools

''' Key for itertools groupby
    Alters flag when sequence changes from annotated from unannotated or vice versa
    Input: genomic sequence
    Output: flag: 0 or 1
'''

class Key(object):
    def __init__(self):
        self.is_nt,self.flag,self.prev = [0,1],[0,1],None
    def __call__(self,e):
        ebool = any(x in self.is_nt for x in e)
        if self.prev:
            prevbool = any(x in self.is_nt for x in self.prev)
        else:
            prevbool = None
        if prevbool != ebool:
            self.flag = self.flag[::-1]
        self.prev = e
        return self.flag[0]


''' fastaN
Return a list of lists separating annotated from unannotated sequence
Input: String
Output: List of lists
Output example: ['ATCG','NNNNNNNNN','AATTTTTTT','N','GCGCGC',...]
'''



def fastaN(fa):
    return [''.join(list(g)) for k,g in itertools.groupby(fa,key=Key())]
