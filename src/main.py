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
ncuts = 50
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
    
    ## If no seq file, it will simulate
    seqs = None
    if len(sys.argv) == 3:
        seqs = aln_reader.read_phylip_cont_file(sys.argv[2])
    else:
        seqs = []
        for i in tree.lvsnms():
            seqs.append(sequence.Sequence(i,""))

    ## Calculate the kernel densities
    x_grid = np.linspace(low, high, ncuts+1)
    for i in seqs:
        #print i.cont_values
        if len(sys.argv) != 3:
            n_basesample = 1000
            x = random.random()
            if x < 0.2:
                i.cont_values = bimodal(0,1,0.5,2,3,2.5,n_basesample)
            else:
                i.cont_values = np.random.rayleigh(1.,n_basesample)
                #i.cont_values = bimodal(1,2,1.5,3,4,3.5,n_basesample)
        density = kde(i.cont_values)
        density.covariance_factor = lambda : .25
        density._compute_covariance()
        kdepdf = density.evaluate(x_grid)
        i.orig_values = i.cont_values
        i.cont_values = kdepdf
    infile.close()
    match_tips_and_cont_values(tree,seqs)
    
    for i in tree.iternodes():
        if len(i.children) != 0:
            i.data['cont_values'] = []

    ## Conduct the analyses
    outfile = open("contanc.tre","w")
    for i in range(ncuts+1):
        curcost = calc_square_change_anc_states(tree,0.00001,i)
        outfile.write(tree.get_newick_repr(True)+"\n")
    outfile.close()

    ## Construct the png plot figures of the tips and internal nodes
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
        plt.figure()
        if len(i.children) > 0:
            i.label = "nd"+str(ndcount)
            ndcount += 1
        else:
            plt.hist(i.data['orig_values'],normed=1,facecolor='g', alpha=0.25,lw=0)
        plt.plot(x_grid,i.data['cont_values'],lw=3,color='blue',alpha=0.5)
        plt.fill_between(x_grid,0,i.data['cont_values'],alpha=0.25, color='blue')
        #plt.show()
        plt.savefig(str(i.label)+'.png')
        plt.close()
        axs[axc][count].plot(i.data['cont_values'])
        axs[axc][count].set_title(i.get_newick_repr(False))
        count += 1
    ts = tree.get_newick_repr(True)+";"
    print ts

    ## Plotting the tree and distribution using ETE
    pl = raw_input("Do you want to plot the tree and results (using ete)? [y/n, default = y]")
    if pl.lower() == 'y' or len(pl) == 0:
        from ete3 import Tree, TreeStyle, TextFace,faces
        from plot_results import mylayout
        t = Tree(ts,format=1)
        # Basic tree style
        ts = TreeStyle()
        ts.show_branch_length=True
        ts.show_leaf_name = True
        ts.layout_fn = mylayout
        # Add two text faces to different columns
        t.show(tree_style=ts)
