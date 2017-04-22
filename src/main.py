import sys,os
import tree_reader,node,sequence
import math
from utils import *
from numpy import *
from anc_state_calculator import *
import aln_reader
import math
from scipy.stats import gaussian_kde as kde
from scipy import stats
import matplotlib.pyplot as plt

#this is the number of cuts
ncuts = 100
#this is the range of the data
low = 0
high = 6
cuts = []
for i in range(ncuts+1):
    cuts.append(low + (i*(float(high-low)/ncuts)))

if __name__ == "__main__":
    if len(sys.argv) != 3 and len(sys.argv) != 2:
        print "python "+sys.argv[0]+" tree [phylipfile]"
        sys.exit(0)
    infile = open(sys.argv[1],"r")
    tree = tree_reader.read_tree_string(infile.readline())
    #if no seq file, it will simulate
    seqs = None
    if len(sys.argv) == 3:
        seqs = aln_reader.read_phylip_cont_file(sys.argv[2])
    else:
        seqs = []
        for i in tree.lvsnms():
            seqs.append(sequence.Sequence(i,""))
    for i in seqs:
        #print i.cont_values
        if len(sys.argv) != 3:
            n_basesample = 1000
            x = random.random()
            if x < 0.2:
                i.cont_values = bimodal(0,1,0.5,2,3,2.5,n_basesample)
            else:
                i.cont_values = np.random.rayleigh(4.,n_basesample)
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
    ndcount = 0
    axc = 0
    for i in tree.iternodes(order="POSTORDER"):
        if count == cn:
            axc += 1
            if axc == 3:
                axc = 0
            count = 0
        if len(i.children) > 0:
            i.label = "nd"+str(ndcount)
            ndcount += 1
        plt.figure()
        plt.plot(i.data['cont_values'])
        plt.savefig(str(i.label)+'.png')
        axs[axc][count].plot(i.data['cont_values'])
        axs[axc][count].set_title(i.get_newick_repr(False))
        count += 1
    #plt.show()
    print tree.get_newick_repr(True)
