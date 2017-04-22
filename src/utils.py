import numpy as np

"""
matches to cont values
"""
def match_tips_and_cont_values(tree,seqs):
    lvs = tree.leaves()
    for i in lvs:
        test = False
        for j in seqs:
            if j.name == i.label:
                i.data['cont_values'] = j.cont_values
                test = True
                break
        if test == False:
            print "can't find "+i.label+" in cont_values"
            return False


def bimodal( low1, high1, mode1, low2, high2, mode2 ,N):
    ns = []
    for i in range(N):
        toss = np.random.choice( (1, 2) )
        if toss == 1:
            ns.append( np.random.triangular( low1, mode1, high1 ) )
        else:
            ns.append( np.random.triangular( low2, mode2, high2 ))
    return ns

def aic(k,l):
    return (2*k)+(2*l)

def aicc(k,l,n):
    return aic(k,l) + ((2*k*(k+1))/(n-k-1))
