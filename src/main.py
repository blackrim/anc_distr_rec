import sys,os
import math
from numpy import *
import math
from scipy.stats import gaussian_kde as kde
from scipy import stats
from datetime import datetime
import argparse
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
plt.style.use('fivethirtyeight')

import tree_reader,node,sequence
from utils import *
from anc_state_calculator import *
import aln_reader



def generate_argparser():
    parser = argparse.ArgumentParser(
        prog="anc_distr_rec",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-t","--treefile",type=open,nargs=1,required=True,
        help=("Input tree in Newick format."))
    parser.add_argument("-d","--datafile",type=open,nargs=1,required=False,
        help=("Datafile in modified Phylip with columns space separated.\
        If not present, then data will be simulated."))
    parser.add_argument("--precut",action="store_true",
        help=("Are the data already in predetermined categories (otherwise, it is expected that the data are independent points)"))
    parser.add_argument("-o","--outdir",type=os.path.abspath,nargs=1,
        required=True,help=("Output directory."))
    parser.add_argument("-c","--ncats",type=int,default=50,
        help=("The number of categories (how much to split it up)."))
    parser.add_argument("--sumtoone",action="store_true",default=True,
        help=("Should we make the distributions sum to one?"))
    parser.add_argument("--showtree",action="store_true",required=False,
        help=("Show the tree (required ete3)."))
    parser.add_argument("--printrate",action="store_true",
        help=("Print the rates in the outdir (requires matplotlib)?"))
    parser.add_argument("--printplots",action="store_true",default=True,
        help=("Print the rates in the outdir (requires matplotlib)?"))
    return parser

def main():
    arguments = sys.argv[1:]
    parser = generate_argparser()
    args = parser.parse_args(arguments)
    #this is the range of the data
    # move this somewhere else soon
    low = 0
    high = 6
    cats = []
    for i in range(args.ncats+1):
        cats.append(low + (i*(float(high-low)/args.ncats)))
    #
    start = datetime.now()
    treeline = args.treefile[0].readline()
    tree = tree_reader.read_tree_string(treeline)
    outd = args.outdir[0]
    if outd[-1] != "/":
        outd += "/"
    print (outd)
    if os.path.isdir(outd) == False:
        os.makedirs(outd)
    ## If no seq file, it will simulate
    seqs = None
    if args.datafile:
        seqs = aln_reader.read_table_cont_file(args.datafile[0])
        if args.precut: 
            ncuts = len(seqs[0].cont_values)
            cuts = list(range(0,args.ncats))
            low = 0
            high = args.ncats
    else:
        seqs = []
        for i in tree.lvsnms():
            seqs.append(sequence.Sequence(i,""))

    ## Calculate the kernel densities
    x_grid = np.linspace(low, high, args.ncats)
    for i in seqs:
        if args.datafile == None:
            n_basesample = 2000
            x = random.random()
            if x < .2:
                i.cont_values = bimodal(0,2,1,4,6,5,n_basesample)
            else:
                i.cont_values = np.random.rayleigh(1.5,n_basesample)
        if args.precut == False or args.datafile == None:
            density = kde(i.cont_values)
            density.covariance_factor = lambda : .25
            density._compute_covariance()
            kdepdf = density.evaluate(x_grid)
            i.orig_values = i.cont_values
            if args.sumtoone:
                i.cont_values,scale = scale_to_one(i.cont_values)
            i.cont_values = kdepdf
        else:
            i.orig_values = i.cont_values
            if args.sumtoone:
                i.cont_values,scale = scale_to_one(i.cont_values)
    match_tips_and_cont_values(tree,seqs)
    

    for i in tree.iternodes():
        if len(i.children) != 0:
            i.data['cont_values'] = []
            i.data['cont_values_se'] = []

    ## Conduct the analyses
    rates = []
    for i in range(args.ncats):
        calc_schluter_anc_states(tree,i)
        estrate = sigsqML(tree)
        rates.append(estrate)

    #plot the estimated rates
    if args.printrate:
        plt.plot(x_grid,rates)
        plt.ylabel('estimated rate')
        plt.xlabel('state space')
        plt.savefig(outd+"rate_estimate.png")
        plt.show()

    ## standard error and scale to one
    for i in tree.iternodes():
        if i.istip == False:
            if args.sumtoone:
                i.data['cont_values'],scale = scale_to_one(i.data['cont_values'])
                i.data['cont_values_se'] = [j*scale for j in i.data['cont_values_se']]
            i.data['cont_values_low'] = [max(j-k,0) for j,k in zip(i.data['cont_values'],i.data['cont_values_se'])]
            i.data['cont_values_high'] = [max(j+k,0) for j,k in zip(i.data['cont_values'],i.data['cont_values_se'])]
    
    end = datetime.now()
    print("Total runtime (H:M:S): "+str(end-start))

    if args.printplots:
        print("creating figures now")
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
                if args.precut == False:
                    plt.hist(i.data['orig_values'],normed=1, histtype='stepfilled',alpha=0.25)
            plt.fill_between(x_grid,0,i.data['cont_values'],alpha=0.05)
            if i.istip == False:
                plt.plot(x_grid,i.data['cont_values_high'],'--',alpha=0.55,)
                plt.plot(x_grid,i.data['cont_values_low'],'--',alpha=0.55,)
            plt.grid(True)
            plt.savefig(outd+str(i.label)+'.png')
            plt.close()

    ts = tree.get_newick_repr(True)+";"

    ## Plotting the tree and distribution using ETE
    if args.showtree and args.printplots:
        from ete3 import Tree, TreeStyle, TextFace,faces,AttrFace
        #from plot_results import mylayout
        def mylayout(node):
            nf = faces.ImgFace(outd+node.name+".png")
            nf.margin_bottom = 20
            nf.margin_right = 20
            nf.margin_left = 20
            faces.add_face_to_node(nf,node,0,position = "branch-top")
            if node.is_leaf():
                ff = AttrFace("name",fsize=30)
                ff.margin_left = 20
                faces.add_face_to_node(ff, node, column=1)
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
        t.render(outd+"contree.svg", w=600, units="mm", tree_style=ts)
        t.render(outd+"contree.pdf", w=6000, units="mm", tree_style=ts)


if __name__ == "__main__":
    main()