import numpy as np
from numpy import *
from scipy.linalg import *

def sigsqML(tree): #tree must already have characters mapped to tips using match_traits_tips()
    n = len(tree.lvsnms())
    vals = [None]*(n-2)
    p = 0
    for i in tree.iternodes():
        i.savelength = i.length
    for i in tree.iternodes(order="POSTORDER"):
        if i.istip == False and i != tree:
            x = [j.data['val'] for j in i.children]
            t = [j.length for j in i.children]
            ui = abs(x[0]-x[1])
            Vi = sum(t)
            vals[p] = (ui,Vi)
            add = (t[0]*t[1])/(t[0]+t[1])
            i.length = i.length + add
            p += 1
        if i == tree:
            t = [j.length for j in i.children]
            Vi = sum(t)
            V0 = (t[0]*t[1])/(t[0]+t[1])
    for i in tree.iternodes():
        i.length = i.savelength
    div = sum([math.pow(i[0],2)/i[1] for i in vals])+(0/V0)
    sig2 = (1./n) * div
    return sig2

"""
calculates square change parsimony for continuous characters
same as likelihood with brownian

NOTE: right now, this isn't using valse (standard error) because we 
are doing distributions
"""
def calc_square_change_anc_states(tree,rate,char, estimaterate = True):
    df = 0
    nodenum = {}
    count = 0
    for i in tree.iternodes(order="postorder"):
        if i.istip:
            i.data['val'] = float(i.data['cont_values'][char])
            i.data['valse'] = float(i.data['cont_values'][char])
        else:
            nodenum[i] = count
            count += 1
            df += 1
            i.data['val'] = 0.
            i.data['valse'] = 0
    df -= 1
    if estimaterate:
        rate = sigsqML(tree)
    #compute the mlest of the root
    fullMcp = zeros((df+1,df+1))
    fullVcp = zeros(df+1)
    count = 0
    for i in tree.iternodes(order="postorder"):
        if i.istip == False:
            nni = nodenum[i]
            for j in i.children:
                tbl = 2./j.length
                fullMcp[nni][nni] += tbl;
                if j.istip:
                    fullVcp[nni] += (j.data['val'] * tbl)
                else:
                    nnj = nodenum[j]
                    fullMcp[nni][nnj] -= tbl;
                    fullMcp[nnj][nni] -= tbl;
                    fullMcp[nnj][nnj] += tbl;
            count += 1
    b = cho_factor(fullMcp)
    #these are the ML estimates for the ancestral states
    mle = cho_solve(b,fullVcp)
    sos = 0
    for i in tree.iternodes(order="postorder"):
        if i.istip == False:
            i.data['val'] = mle[nodenum[i]]
            i.data['cont_values'].append(mle[nodenum[i]])
            #print i.data['val']
            #i.label = str(mle[nodenum[i]])
            for j in i.children:
                temp = (i.data['val'] - j.data['val'])
                sos += temp*temp / j.length
    #print "Square Length: ",sos
    #calcSE
    for i in tree.iternodes(order="postorder"):
        if i.istip == False:
            qpq = fullMcp[nodenum[i]][nodenum[i]]
            tm1 = np.delete(fullMcp,(nodenum[i]),axis=0)
            tm = np.delete(tm1,(nodenum[i]),axis=1)
            b = cho_factor(tm)
            sol = cho_solve(b,tm1[:,nodenum[i]])
            tempse = qpq - np.inner(tm1[:,nodenum[i]],sol)
            i.data['valse'] = math.sqrt(2*sos/(df*tempse))
            i.data['cont_values_se'].append(i.data['valse'])
    return rate

