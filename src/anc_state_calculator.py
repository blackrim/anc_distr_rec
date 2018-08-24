import numpy as np
import tree_utils
from numpy import *
from scipy.linalg import *

def var_cov_matrix(tree):
    lvs = tree.leaves()
    ll = len(lvs)
    zm = zeros((ll,ll))
    for i in range(ll):
        for j in range(ll):
            curnode = None
            if i == j:
                curnode = lvs[i]
            elif j > i:
                curnode = tree_utils.get_mrca([lvs[i],lvs[j]],tree)
            else:
                continue
            count = 0
            while curnode != tree:
                count += curnode.length
                curnode = curnode.parent
            zm[i][j] = count
            zm[j][i] = count
    return zm

def sigsqML(tree):
    lvs = tree.leaves()
    ahat = tree.data['val']
    N = len(lvs)
    EX = np.zeros((N,1))
    for i in range(N):
        EX[i][0] = float(tree.data['val'])
    X = np.zeros((N,1))
    for i in range(N):
        X[i][0] = float(lvs[i].data['val'])
    ones = np.ones((N,1))
    vcv = var_cov_matrix(tree)
    return (((X-EX).T.dot(inv(vcv)).dot(X-EX))/(N))[0]

def calc_schluter_anc_states(tree,char):
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
    # NOTE: these two cholesky steps are the slowest bits
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
            #print i.get_newick_repr(False),mle[nodenum[i]]
            for j in i.children:
                temp = (i.data['val'] - j.data['val'])
                sos += temp*temp / j.length
    #calcSE
    # need to change this to the rohlf 2001 standard errors
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
            #print i.data['valse']

