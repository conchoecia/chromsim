"""
This file contains various utility functions for the program, like plotting and logging.
"""

import os
import matplotlib.pyplot as plt
from scipy import stats
import numpy as np
import collections.abc as abc
import argparse as ap
import yaml
import glob
from datetime import datetime as dt
import pandas as pd
from chromosome import Chrom

csv_separator='\t'

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
    
def init_plot_style_settings():
    """
    set the global values for text sizes, colors, etc. to be used in plots
    """

    # figure styling
    plt.style.use('bmh')

    # set the colors, font sizes, etc. as global variables
    global setlw
    setlw = 0.2
    global seta
    seta  = 0.1
    global fc
    fc='white'
    global box_alpha
    box_alpha=0.5
    global ec
    ec='gray'
    global bbox
    bbox=dict(facecolor=fc, alpha=box_alpha, boxstyle='round, pad=0.5', edgecolor=ec)
    global plot_title_size
    plot_title_size=50
    global subplot_title_size
    subplot_title_size=30
    global text_size
    text_size=20    

def set_up_dot_fig(chrom, title):
    """
    set up all the axes needed for the dotplots in a grid
    """
    
    figsize=(40, 40)
    fig=plt.figure(figsize=figsize)
    gs=fig.add_gridspec(2, 2)
    ax0=fig.add_subplot(gs[0, 0])
    ax1=fig.add_subplot(gs[0, 1])
    ax2=fig.add_subplot(gs[1, 0])
    ax3=fig.add_subplot(gs[1, 1])

    subplot0_title=r"$\tau_0$ vs. $\tau_{100\%}$"
    subplot1_title=r"$\tau_0$ vs. $\tau_{50\%}$"
    subplot2_title=r"$\tau_0$ vs. $\tau_S$"
    subplot3_title=r"$\tau_S$ vs. $\tau_{50\%}$"

    ax0.set_xlabel(r"$\tau_0$", fontsize=text_size)
    ax0.set_ylabel(r"$\tau_{100\%}$", fontsize=text_size)
    ax1.set_xlabel(r"$\tau_0$", fontsize=text_size)
    ax1.set_ylabel(r"$\tau_{50\%}$", fontsize=text_size)
    ax2.set_xlabel(r"$\tau_0$", fontsize=text_size)
    ax2.set_ylabel(r"$\tau_S$", fontsize=text_size)
    ax3.set_xlabel(r"$\tau_S$", fontsize=text_size)
    ax3.set_ylabel(r"$\tau_{50\%}$", fontsize=text_size)

    fig.suptitle(title, fontsize=plot_title_size)
    ax0.set_title(subplot0_title, fontsize=subplot_title_size)
    ax1.set_title(subplot1_title, fontsize=subplot_title_size)
    ax2.set_title(subplot2_title, fontsize=subplot_title_size)
    ax3.set_title(subplot3_title, fontsize=subplot_title_size)

    return fig, ax0, ax1, ax2, ax3
    
def set_up_trace_fig(chrom, title):
    """
    set up all the axes needed for the traces on a grid
    """

    figsize=(20, 25)
    fig=plt.figure(figsize=figsize)
    gs=fig.add_gridspec(3, 2, width_ratios=[3, 1])
    ax0=fig.add_subplot(gs[0, :]) # add subplot for gene interaction traces spanning the top half
    ax1=fig.add_subplot(gs[2, 0]) # add subplot for m values
    ax2=fig.add_subplot(gs[2, 1], sharey=ax1) # add subplot for m value normal distribution
    ax3=fig.add_subplot(gs[1, 0], sharex=ax1) # add subplot for interactions in the range of the m plot
    ax0.set_xlim(0, chrom.cycle)
    ax1.set_xlim([0, chrom.tS*2])
    
    global A_alpha
    A_alpha   = max(1/chrom.genesA, 0.1)
    global B_alpha
    B_alpha   = max(1/chrom.genesB, 0.1)
    global all_alpha
    all_alpha = max(1/(chrom.genesA + chrom.genesB), 0.05)
    
    # set plot titles and axis labels
    subplot0_title=r"number of unique interactions after $n$ inversion cycles"
    subplot1_title=r"$m$ after $n$ inversion cycles"

    ax0.set_xlabel("inversion cycle", fontsize=text_size)
    ax0.set_ylabel("unique interactions", fontsize=text_size)
    ax1.set_xlabel("inversion cycle", fontsize=text_size)
    ax1.set_ylabel(r"$m$", fontsize=text_size)
    ax3.set_xlabel("inversion cycle", fontsize=text_size)
    ax3.set_ylabel("unique interactions", fontsize=text_size)

    fig.suptitle(title, fontsize=plot_title_size)
    ax0.set_title(subplot0_title, fontsize=subplot_title_size)
    ax1.set_title(subplot1_title, fontsize=subplot_title_size)
    
    return fig, ax0, ax1, ax2, ax3

def plot_dotplot(chrom, ax, x, y):
    """
    plot a dotplot between two gene lists
    """

    # get the indices of all genes in x in both lists
    # segregated into the two gene groups to make them different color dots
    xA_indices=[]
    yA_indices=[]
    xB_indices=[]
    yB_indices=[]
    for gene in x:
        if gene.startswith('A'):
            xA_indices.append(x.index(gene))
            yA_indices.append(y.index(gene))
        elif gene.startswith('B'):
            xB_indices.append(x.index(gene))
            yB_indices.append(y.index(gene))

    ax.scatter(xA_indices, yA_indices, lw=setlw*2)
    ax.scatter(xB_indices, yB_indices, lw=setlw*2)

def plot_trace(chrom, trace, ax, lim, color='blue', alpha=0.5):
    """
    plot trace dictionary on the given axes in the interval [0, lim]
    """

    # plot all the individual traces
    # at the same time collect them into a single combined trace
    all=[]
    for k in trace:
        x, y= [], []
        for tup in trace[k]:
            if tup[0] > lim:
                break
            x.append(tup[0])
            y.append(tup[1])
            all.append(tup)
        ax.plot(x, y, color=color, lw = setlw, alpha=alpha)
    all=sorted(all)

    # calculate the moving average of the combined trace
    df = pd.DataFrame(all, columns =['cycle', 'interactions'])
    new_all=df.rolling(len(all)//10).mean()

    # add the last elements of the original list to make sure the line goes all the way to the end
    x=new_all['cycle'].tolist()
    x.append(all[-1][0])
    y=new_all['interactions'].tolist()
    y.append(all[-1][1])

    # plot the combined trace
    ax.plot(x, y, color=color, lw = setlw*10, alpha=0.75)

def plot_t50(chrom, ax):
    """
    plot a horizontal line at x=t50 on the given axis
    """
    
    t50_text=r"$\tau_{{50\%}}={t50}$" "\n" "$({perc:.2f}\%\ of\ cycles)$".format(t50=chrom.t50, perc=chrom.t50/chrom.cycle*100)
    ax.text(x=chrom.t50, y=1, ha='left', va='center', s=t50_text, bbox=bbox, fontsize=text_size)
    ax.axvline(x=chrom.t50, lw=setlw*5, color='black')

def plot_m(chrom, ax_m, ax_norm):
    """
    plot the m value curve and the normal distribution (on a separate plot ax_norm)

         ax_m    ax_norm
    +-----------+ +---+
    |   /-------| |\  |
    |  /        | | | |
    | /         | |/  |
    |/          | |   |
    +-----------+ +---+
    """

    cycles=[x for x in chrom.trace_m.keys()]
    m_values=[y for y in chrom.trace_m.values()]
    pdf_space=np.linspace(0, max(m_values), 100)
    normpdf=stats.norm.pdf(pdf_space, chrom.m_mu, chrom.m_sigma)
    normpdf/=max(normpdf)
    upper_bound=chrom.m_mu+1.96*chrom.m_sigma
    lower_bound=chrom.m_mu-1.96*chrom.m_sigma
        
    ax_m.plot(cycles[:chrom.tS*2+1], m_values[:chrom.tS*2+1], lw=setlw*2, color='blue', label=r"$m$")
    
    # plot 95 percentile of m value normal distribution
    crossed_text="first 95 percentile value:\n{cross} cycles\n({perc:.2f}% of t50)".format(cross=chrom.tS, perc=chrom.tS/chrom.t50*100)
    norm_label=r"normal distribution of $m$" "\n" "(excluding the first {perc}% of cycles)".format(perc=25)
    ax_m.axhline(y=upper_bound, color='red', lw=setlw*5, ls=':') # plot upper bound of the 95 percentile
    ax_m.axhline(y=lower_bound, color='red', lw=setlw*5, ls=':') # plot lower bound of the 95 percentile
    ax_m.text(y=upper_bound, x=chrom.tS*2, ha='right', va='bottom', s=r"$\mu+1.96\cdot\sigma$", bbox=bbox, fontsize=text_size)
    ax_m.text(y=lower_bound, x=chrom.tS*2, ha='right', va='top', s=r"$\mu-1.96\cdot\sigma$", bbox=bbox, fontsize=text_size)
    ax_m.fill_between(cycles, lower_bound, upper_bound, color='red', alpha=0.1) # shade area between the bounds of the 95 percentile
    ax_m.axvline(x=chrom.tS, lw=setlw*5, color='black') # plot the x value where the m value first enters the 95 percentile
    ax_m.text(x=chrom.tS, y=0.4, ha='left', va='center', s=crossed_text, bbox=bbox, fontsize=text_size)
    
    # plot normal distribution next to m plot
    ax_norm.plot(normpdf, pdf_space, lw=setlw*5, color='red', label=norm_label) # plot the normal distribution of the m values along the y axis

def save_fig(output_dir, output_name):
    """
    save the figure of the current matplotlib plot to a specified output location as .png and .pdf
    """
    
    # create the diagram directory if it does not exist yet
    diagram_dir=output_dir+'diagrams/'
    
    if not os.path.exists(diagram_dir):
        os.makedirs(diagram_dir)
        if not os.path.exists(diagram_dir+'/pdf'):
            os.makedirs(diagram_dir+'/pdf')
            if not os.path.exists(diagram_dir+'/png'):
                os.makedirs(diagram_dir+'/png')
                
    # save this as a pdf and png
    plt.savefig(diagram_dir+'pdf/'+output_name+'.pdf')
    plt.savefig(diagram_dir+'png/'+output_name+'.png')
    
    #if not yaml:
    #    return
    #
    ## save the trace as a yaml file
    #if not os.path.exists(diagram_dir+'/yaml'):
    #    os.makedirs(diagram_dir+'/yaml')
    #    with open(diagram_dir+'yaml/'+output_name+'.yaml', 'w') as f:
    #        yaml.dump(chrom.trace, f)

def plot_results(chrom, output_dir):
    """
    plot the results of a simulated chromosome
    """    

    init_plot_style_settings()

    output_name=get_output_name(chrom)
    plot_title=get_plot_title(chrom)
    
    ## trace figure

    # calculate the upper bounds up to which should be plotted
    m_lim=(chrom.tS*2)
    lim=(chrom.cycle)

    # set up the fig
    fig, ax0, ax1, ax2, ax3=set_up_trace_fig(chrom, plot_title)

    # plot traces
    plot_trace(chrom, chrom.trace, ax0, lim, 'black', all_alpha)
    plot_trace(chrom, chrom.trace, ax3, m_lim, 'black', all_alpha)
    plot_trace(chrom, chrom.trace_AtoB, ax0, lim, 'red', A_alpha)
    plot_trace(chrom, chrom.trace_AtoB, ax3, m_lim, 'red', A_alpha)
    plot_trace(chrom, chrom.trace_BtoA, ax0, lim, 'blue', B_alpha)
    plot_trace(chrom, chrom.trace_BtoA, ax3, m_lim, 'blue', B_alpha)  

    plot_t50(chrom, ax0)
    plot_m(chrom, ax1, ax2)

    # set legends
    ax1.legend(facecolor=fc, framealpha=box_alpha, edgecolor=ec)
    ax2.legend(facecolor=fc, framealpha=box_alpha, edgecolor=ec)

    # save the trace figure
    save_fig(output_dir, output_name+"_trace")

    ## dotplot figure

    # get the gene lists
    original_genes=chrom.original_gene_list
    final_genes=chrom.gene_list
    half_time_genes=chrom.t50_gene_list

    # rerun the chromosome until tS to get the list at that point
    ts_chrom=Chrom(chrom.length, chrom.genesA, chrom.genesB, chrom.level_of_convergence, chrom.window_size, chrom.translocations_per_cycle, chrom.inversion_cuts)
    ts_chrom.run(n=chrom.tS, show_output=False)
    entropy_genes=ts_chrom.gene_list

    # set up the fig
    fig, ax0, ax1, ax2, ax3=set_up_dot_fig(chrom, plot_title)

    # plot dotplots
    plot_dotplot(chrom, ax0, original_genes, final_genes)
    if chrom.t50 >= 0:
        plot_dotplot(chrom, ax1, original_genes, half_time_genes)
    if chrom.tS >= 0:
        plot_dotplot(chrom, ax2, original_genes, entropy_genes)
    if chrom.tS >= 0 and chrom.t50 >= 0:
        plot_dotplot(chrom, ax3, entropy_genes, half_time_genes)

    # save the dotplot figure
    save_fig(output_dir, output_name+"_dotplots")

def log(chrom, output_dir, elapsed='N/A'):
    """
    log the state of a chromosome

    TODO: change this to a better format
    """

    # write the state to a log file according to the input
    log_dir=output_dir+'log/'
    log_file=get_output_name(chrom)+'.csv'

    if not os.path.exists(log_dir):
        os.makedirs(log_dir)

    new_log_file=not os.path.exists(log_dir+log_file)
    mode='w' if new_log_file else 'a'
    
    log_header='timestamp;|A|;|B|;cycles;t50;AB_convergence;first_95_m;m_sigma;m_mu;Delta_t\n'
    log_format_string='{ts};{A};{B};{c};{t50};{ABconv};{m95};{sig:.3f};{mu:.3f};{dt}\n'
    log_line=log_format_string.format(ts=str(dt.now()), A=chrom.genesA, B=chrom.genesB, c= chrom.cycle, t50=chrom.t50, ABconv=chrom.t100, m95=chrom.tS, sig=chrom.m_sigma, mu=chrom.m_mu, dt=elapsed)

    with open(log_dir+log_file, mode) as f:
        if new_log_file:
            f.write(log_header)
        f.write(log_line)

def metalog(output_dir):
    """
    TODO: better format
    """
    
    # write some information about the input paramters into a separate log file
    log_dir=output_dir+'log/'
    metalog_file='metalog.csv'

    metalog_header='|A|;|B|;converging;window_size;level_of_convergence;runs;average_t50;average_t100;average_ts;average_runtime\n'
    metalog_format_string='{A};{B};{conv};{wsize};{loc};'
    metainfo_format_string='{runs};{avt50:.2f};{avt100:.2f};{avts:.2f};{avrt}\n'

    logfiles=glob.glob(log_dir+'inversion_sim*')
    with open(log_dir+metalog_file, 'w') as mf:
        mf.write(metalog_header)

        for logfile in logfiles:
            s=logfile.split('/')[-1].split('_')
            sizes=s[2].split('-')
            Asize=int(sizes[0][1:])
            Bsize=int(sizes[1][1:])
            loc=int(s[3][1:])
            wsize=int(s[4][1:])
            converging=(len(s) == 6)
            runs=[]

            with open(logfile) as lf:
                runs=[run.rstrip().split(csv_separator) for run in lf.readlines()]

            header=runs[0]
            t50s=[float(run[header.index('t50')]) for run in runs[1:]]
            t100s=[float(run[header.index('cycles')]) for run in runs[1:]]
            tss=[float(run[header.index('first_95_m')]) for run in runs[1:]]
            rts=[runtime_to_int(run[header.index('Delta_t')]) for run in runs[1:]]
            
            metalog_line=metalog_format_string.format(A=Asize, B=Bsize, conv=converging, wsize=wsize, loc=loc)        
            metainfo_string=metainfo_format_string.format(runs=len(runs), avt50=np.average(t50s), avt100=np.average(t100s), avts=np.average(tss), avrt=int_to_runtime(np.average(rts)))
            mf.write(metalog_line+metainfo_string)

def get_output_name(chrom):
    """
    return a string to be used as the base for saving logs

    TODO: will be redundant with the new log file format
    """
    
    return 'inversion_sim_A{Asize}-B{Bsize}_l{loc}_w{wsize}'.format(Asize=chrom.genesA, Bsize=chrom.genesB, loc=chrom.level_of_convergence, wsize=chrom.window_size)

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

