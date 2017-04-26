import sys,os
import math
from numpy import *
import math
from scipy.stats import gaussian_kde as kde
from scipy import stats
from datetime import datetime
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
plt.style.use('fivethirtyeight')

import tree_reader,node,sequence
from utils import *
from anc_state_calculator import *
import aln_reader
import conf

#this is the number of cuts
ncuts = 25
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
    start = datetime.now()
    infile = open(sys.argv[1],"r")
    tree = tree_reader.read_tree_string(infile.readline())
    
    if os.path.isdir(conf.outdir) == False:
        os.makedirs(conf.outdir)

    ## If no seq file, it will simulate
    seqs = None
    if len(sys.argv) == 3:
        seqs = aln_reader.read_table_cont_file(sys.argv[2])
        if conf.precut: 
            ncuts = len(seqs[0].cont_values)
            cuts = range(0,ncuts)
            low = 0
            high = ncuts
    else:
        seqs = []
        for i in tree.lvsnms():
            seqs.append(sequence.Sequence(i,""))

    ## Calculate the kernel densities
    x_grid = np.linspace(low, high, ncuts)
    for i in seqs:
        if len(sys.argv) != 3:
            n_basesample = 2000
            x = random.random()
            if x < 0.2:
                i.cont_values = bimodal(0,2,1,4,6,5,n_basesample)
            else:
                i.cont_values = np.random.rayleigh(1.5,n_basesample)
        if conf.precut == False:
            density = kde(i.cont_values)
            density.covariance_factor = lambda : .25
            density._compute_covariance()
            kdepdf = density.evaluate(x_grid)
            i.orig_values = i.cont_values
            if conf.sumtoone:
                i.cont_values,scale = scale_to_one(i.cont_values)
            i.cont_values = kdepdf
        else:
            i.orig_values = i.cont_values
            if conf.sumtoone:
                i.cont_values,scale = scale_to_one(i.cont_values)
    infile.close()
    match_tips_and_cont_values(tree,seqs)
    
    for i in tree.iternodes():
        if len(i.children) != 0:
            i.data['cont_values'] = []
            i.data['cont_values_se'] = []

    ## Conduct the analyses
    rates = []
    for i in range(ncuts):
        estrate = calc_square_change_anc_states(tree,0.00001,i)
        print estrate
        rates.append(estrate)

    #plot the estimated rates
    plt.plot(x_grid,rates)
    plt.ylabel('estimated rate')
    plt.xlabel('state space')
    plt.savefig(conf.outdir+"rate_estimate.png")
    plt.show()

    ## standard error and scale to one
    for i in tree.iternodes():
        if i.istip == False:
            if conf.sumtoone:
                i.data['cont_values'],scale = scale_to_one(i.data['cont_values'])
                i.data['cont_values_se'] = [j*scale for j in i.data['cont_values_se']]
            i.data['cont_values_low'] = [max(j-k,0) for j,k in zip(i.data['cont_values'],i.data['cont_values_se'])]
            i.data['cont_values_high'] = [max(j+k,0) for j,k in zip(i.data['cont_values'],i.data['cont_values_se'])]
    
    end = datetime.now()
    print "Total runtime (H:M:S): "+str(end-start)
    print "creating figures now"
    ## Construct the png plot figures of the tips and internal nodes
    ndcount = 0
    for i in tree.iternodes(order="POSTORDER"):
        plt.figure(figsize=(6, 4))
        plt.plot(x_grid,i.data['cont_values'])
        #should output the i.data['cont_values'] for each node here
        if len(i.children) > 0:
            i.label = "nd"+str(ndcount)
            ndcount += 1
        else:
            if conf.precut == False:
                plt.hist(i.data['orig_values'],normed=1, histtype='stepfilled',alpha=0.25)
        plt.fill_between(x_grid,0,i.data['cont_values'],alpha=0.05)
        if i.istip == False:
            plt.plot(x_grid,i.data['cont_values_high'],'--',alpha=0.55,)
            plt.plot(x_grid,i.data['cont_values_low'],'--',alpha=0.55,)
        plt.grid(True)
        plt.savefig(conf.outdir+str(i.label)+'.png')
        plt.close()
    ts = tree.get_newick_repr(True)+";"

    ## Plotting the tree and distribution using ETE
    if conf.showtree:
        from ete3 import Tree, TreeStyle, TextFace,faces
        from plot_results import mylayout
        t = Tree(ts,format=1)
        # Basic tree style
        ts = TreeStyle()
        ts.show_branch_length=False
        ts.show_leaf_name = False
        ts.scale =  900
        ts.branch_vertical_margin = 1
        ts.layout_fn = mylayout
        # Add two text faces to different columns
        #t.show(tree_style=ts)
        t.render("contree.svg", w=600, units="mm", tree_style=ts)
        t.render("contree.pdf", w=6000, units="mm", tree_style=ts)
