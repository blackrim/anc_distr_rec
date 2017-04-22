import sys,os
import tree_reader,node
import math
import numpy as np
from numpy import *
from anc_state_calculator import *
import aln_reader
import math
from scipy.stats import gaussian_kde as kde
from scipy import stats
import matplotlib.pyplot as plt

ncuts = 100
low = 0
high = 6
cuts = []
for i in range(ncuts+1):
    cuts.append(low + (i*(float(high-low)/ncuts)))

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


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print "python "+sys.argv[0]+" tree phy"
        sys.exit(0)
    infile = open(sys.argv[1],"r")
    tree = tree_reader.read_tree_string(infile.readline())
    seqs = aln_reader.read_phylip_cont_file(sys.argv[2])
    for i in seqs:
        #print i.cont_values
        n_basesample = 1000
        x = random.random()
        if x < 0.2:
            i.cont_values = bimodal(0,1,0.5,2,3,2.5,n_basesample)
        else:
            i.cont_values = np.random.rayleigh(0.5,n_basesample)
        density = kde(i.cont_values)
        x_grid = np.linspace(low, high, ncuts+1)
        kdepdf = density.evaluate(x_grid)
        #plt.figure()
        ##density.covariance_factor = lambda : .25
        ##density._compute_covariance()
        #plt.hist(i.cont_values,normed=1)
        #plt.plot(x_grid,kdepdf)
        #plt.xlabel(i.name)
        #plt.show()
        #print kdepdf
        i.cont_values = kdepdf
    infile.close()
    match_tips_and_cont_values(tree,seqs)
    #calc the parsimony
    
    for i in tree.iternodes():
        if len(i.children) != 0:
            i.data['cont_values'] = []

    outfile = open("contanc.tre","w")
    for i in range(ncuts+1):
        curcost = calc_square_change_anc_states(tree,0.00001,i)
        outfile.write(tree.get_newick_repr(True)+"\n")
    outfile.close()

    totaln = len(list(tree.iternodes()))
    cn = int(round((totaln)/3+1))
    f, axs = plt.subplots(3,cn, sharex=False, sharey=False)
    count = 0
    axc = 0
    for i in tree.iternodes(order="POSTORDER"):
        if count == cn:
            axc += 1
            if axc == 3:
                axc = 0
            count = 0
        axs[axc][count].plot(i.data['cont_values'])
        axs[axc][count].set_title(i.get_newick_repr(False))
        count += 1
    plt.show()
