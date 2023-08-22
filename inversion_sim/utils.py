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
        
def save(chrom, outdir="./"):
    """
    save the chromosome in a re-runnable file format
    """
        
    header_format_string="""# timestamp {ts}
## params
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

    fname=str(chrom.timestamp)+FILE_ENDING

    header=header_format_string.format(ts=chrom.timestamp, Asize=chrom.Asize, Bsize=chrom.Bsize, ws=chrom.window_size, loc=chrom.level_of_convergence, t100=chrom.t100, t50=chrom.t50, tS=chrom.tS, mu=chrom.m_mu, sigma=chrom.m_sigma, AB_convergence=chrom.AB_convergence)

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
            
        # extract the parameters and results from the header
        params=[line.split(' ')[2] for line in head[2:6]]
        results=[line.split(' ')[2] for line in head[7:13]]

        # parse the cuts to a list of tuples
        cuts=[]
        line=f.readline()
        while line != '':
            ints=list(map(int, line.rstrip().split('\t')))
            cuts.append((ints[0], ints[1]))
            line=f.readline()

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
        ts=head[0]

        return chrom, results_dict, ts

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

def create_parser():
    """
    creates an argparse parser for the CLI arguments
    """
    parser=ap.ArgumentParser(prog="inversion_sim", description="This program simulates inversion events of a chromosome made up of A and B genes")

    # simulation arguments
    parser.add_argument('-S', '--simulate', action='store_true', help="simulate a chromosome with the parameters -a, -b, -c, -l, -w, and -t")
    parser.add_argument('-a', '--asize', type=int, default=-1, help="integer value for the number of genes in group A")
    parser.add_argument('-b', '--bsize', type=int, default=-1, help="integer value for the number of genes in group B")
    parser.add_argument('-c', '--cycle-number', type=int, default=-1, help="integer value for the number of cycles to run (optional)")
    parser.add_argument('-l', '--level-of-convergence', type=float, metavar='LOC', choices=FloatRange(0, 10), default=1, help="fraction of possible gene interactions to wait for if converging (optional)")
    parser.add_argument('-w', '--window-size', type=int, default=1, help="the size of the window to the left and right of each gene to count as interaction after each cycle (optional)")
    parser.add_argument('-t', '--translocations-per-cycle', type=int, default=0, help="integer value for the number of translocations to be done in addition to inversion each cycle (optional)")

    # plotting arguments
    parser.add_argument('-P', '--plot', action='store_true', help="plot a chromosome specified by -s")
    parser.add_argument('-s', '--source', type=str, help="the .inv file to load")
    parser.add_argument('-G', '--gif', action='store_true', help="create an animated GIF of the simulation process")
    
    # other
    parser.add_argument('-o', '--output-dir', default=os.getcwd(), help="directory in which to store the output of the program (optional)")
    
    return parser

