import numpy as np

"""
matches to cont values
"""
def match_tips_and_cont_values(tree,seqs):
    ln = {}
    sq = {}
    for i in tree.leaves():
        ln[i.label] = i
    for i in seqs:
        sq[i.name] = i
    for i in ln:
        try:
            ln[i].data['cont_values'] = sq[i].cont_values
            ln[i].data['orig_values'] = sq[i].orig_values
        #test = False
        #for j in seqs:
        #    if j.name == i.label:
        #        i.data['cont_values'] = j.cont_values
        #        i.data['orig_values'] = j.orig_values
        #        test = True
        #        break
        #if test == False:
        except:
            print("can't find "+i.label+" in cont_values")
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

def scale_to_one(li):
    sc = float(1)/sum(li)
    nl = []
    for i in li:
        nl.append(i*sc)
    return nl,sc

def aic(k,l):
    return (2*k)+(2*l)

def aicc(k,l,n):
    return aic(k,l) + ((2*k*(k+1))/(n-k-1))
