"""
This file contains various utility functions for the program, like plotting and logging.
"""

import os
import collections.abc as abc
import argparse as ap
import glob
from datetime import datetime as dt
from chromosome import Chrom

FILE_ENDING=".inv"

def get_outname(chrom):
    """
    return an output name for file related to the chromosome
    """
    
    return str(chrom.timestamp)
        
def save(chrom, outdir="./"):
    """
    save the chromosome in a re-runnable file format
    """
        
    header_format_string="""## params
# Asize {Asize}
# Bsize {Bsize}
# window_size {ws}
# level_of_convergence {loc:.2f}
## results
# t100 {t100}
# t50 {t50}
# tS {tS}
# mu {mu:.2f}
# sigma {sigma:.2f}
# AB_convergence {AB_convergence}

"""

    fname=get_outname(chrom)+FILE_ENDING

    header=header_format_string.format(Asize=chrom.Asize, Bsize=chrom.Bsize, ws=chrom.window_size, loc=chrom.level_of_convergence, t100=chrom.t100, t50=chrom.t50, tS=chrom.tS, mu=chrom.m_mu, sigma=chrom.m_sigma, AB_convergence=chrom.AB_convergence)

    with open(outdir+fname, 'w') as f:
        f.write(header)
        for cut in chrom.inversion_cuts:
            f.write('{}\t{}\n'.format(cut[0], cut[1]))

def parse_inv_file(file):
    """
    parse a chromosome from a .inv file
    """

    with open(file, 'r') as f:
        # get header of the files
        head=[]
        line=f.readline()
        while line != '\n':
            head.append(line.rstrip())
            line=f.readline()
            
        print("head")
        print(head)
        # extract the parameters and results from the header
        params=[line.split(' ')[2] for line in head[1:5]]
        results=[line.split(' ')[2] for line in head[6:12]]

        # parse the cuts to a list of tuples
        cuts=[]
        while line !='\n':
            line=f.readline()
            ints=list(map(int, line.rstrip().split('\t')))
            cuts.append((ints[0], ints[1]))

        # parse the parameters and results
        Asize=int(params[0])
        Bsize=int(params[1])
        ws=int(params[2])
        loc=float(params[3])
        
        t100=int(results[0])
        t50=int(results[1])
        tS=int(results[2])
        mu=float(results[3])
        sigma=float(results[4])
        AB_conv=int(results[5])

        # set up the return values
        chrom=Chrom(Asize=Asize, Bsize=Bsize, level_of_convergence=loc, window_size=ws, inversion_cuts=cuts)
        results_dict={'t100': t100,
                      't50': t50,
                      'tS': tS,
                      'mu': mu,
                      'sigma': sigma,
                      'AB_convergence': AB_conv}

        return chrom, results_dict

class FloatRange(abc.Container):

    def __init__(self, lower, upper, step=0.01):
        self.lower=lower
        self.upper=upper
        # step is for the iterator
        self.step=step

    def __contains__(self, x):
        return self.lower <= x and self.upper >= x

    # for some reason this is needed to have a help message with argparse
    def __iter__(self):
        self.n=0
        return self
    
    def __next__(self):
        if self.n>=1:
            raise StopIteration
        self.n+=self.step
        return self.n-self.step

def get_plot_title(chrom):
    """
    return a string to be used as title for the plots of this chromosome

    TODO: make a bit better (idk how yet)
    """
    return r"$|A|={Asize}, |B|={Bsize}, w={wsize}$".format(Asize=chrom.genesA, Bsize=chrom.genesB, wsize=chrom.window_size)

def create_parser():
    """
    creates an argparse parser for the CLI arguments
    """
    parser=ap.ArgumentParser(prog="inversion_sim", description="This program simulates inversion events of a chromosome made up of A and B genes")

    # positional arguments
    parser.add_argument('Asize', type=int, nargs='?', help="integer value for the number of genes in group A")
    parser.add_argument('Bsize', type=int, nargs='?', help="integer value for the number of genes in group B")

    # named general arguments
    parser.add_argument('-o', '--output-dir', default=os.getcwd(), help="directory in which to store the output of the program (default: './')")

    # named simulation arguments
    parser.add_argument('-c', '--cycle-number', type=int, default=-1, help="integer value for the number of cycles to run (alternative to --converge)")
    parser.add_argument('-l', '--level-of-convergence', type=float, metavar='LOC', choices=FloatRange(0, 10), default=1, help="fraction of possible gene interactions to wait for if converging (default: 1)")
    parser.add_argument('-w', '--window-size', type=int, default=1, help="the size of the window to the left and right of each gene to count as interaction after each cycle (default: 1)")
    parser.add_argument('-t', '--translocations-per-cycle', type=int, default=0, help="integer value for the number of translocations to be done in addition to inversion each cycle (default: 0)")

    # flag simulation arguments
    parser.add_argument('-P', '--plot-curves', action='store_true', default=False, help="tell the simulation to plot the curves in the end (default: False)")

    # flag simulation-free arguments
    parser.add_argument('-M', '--metalog', action='store_true', help="collect averages from log files and store them in metalog.csv (not running a simulation)")
    parser.add_argument('-T', '--plot-average-t50', action='store_true', help="plot the average t50 of all previous simulations in this directory (not running a simulation)")
    
    return parser

