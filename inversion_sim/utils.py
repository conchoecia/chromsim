import os
import matplotlib.pyplot as plt
from scipy import stats
import numpy as np
import collections.abc as abc
import argparse as ap
import yaml
from datetime import datetime as dt

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

def convert_dt(dt):
    """
    converts a sting of the form 'mm:ss' to seconds
    """
    mins, secs=dt.split(':')
    return 60*int(mins)+int(secs)

def read_log_file(path, filename, cols=[]):
    """
    read a log file and return the desired columns (all if col is [])
    """
    path+='log/'
    input_file=filename+'.csv'

    if not os.path.exists(path+input_file):
        print(path+input_file)
        print('log file does not exist')
        return

    data=None
    
    with open(path+input_file) as f:
        lines=[line.rstrip().split(';') for line in f.readlines()]
        header=lines[0]
        if cols == []:
            cols=header
        indices=[header.index(col) for col in cols]
        data=[[line[i] for i in indices] for line in lines[1:]]

    return data

#def read_metalog_file(path, cols=[]):
#    """
#    read the metalog file and return the desired columns (all if col is [])
#    """
#    path+='log/'
#    input_file='metalog.csv'
#
#    if not os.path.exists(path+input_file):
#        print(path+input_file)
#        print('metalog file does not exist')
#        return

    

def plot_runtime(path, filename):
    """
    plot the runtime as a function of chromosome size from a log file
    probably obsolete now, will keep just in case
    """
    # might still be useful, but I got to change the code since the log files are named differently now
    output_file=filename+'.png'
    raw_data=read_log_file(path, filename, ['|A|', '|B|', 'Delta_t'])
    data=[[line[0]+line[1], line[2]] for line in raw_data]
    
    plt.plot([x[0] for x in data], [x[1]/60 for x in data], 'o')
    plt.xlabel("chromosome size (|A|+|B|)")
    plt.ylabel("runtime (minutes)")
    plt.savefig(path+output_file)

def set_up_fig(chrom):
    """
    create the basic structure of the matplotlib figure
    returns the figure itself and the axes that are part of it
    """
    import matplotlib.pyplot as plt

    figsize=(20, 25)
    fig=plt.figure(figsize=figsize)
    gs=fig.add_gridspec(3, 2, width_ratios=[3, 1])
    ax0=fig.add_subplot(gs[0, :]) # add subplot for gene interaction traces spanning the top half
    ax1=fig.add_subplot(gs[2, 0]) # add subplot for m values
    ax2=fig.add_subplot(gs[2, 1], sharey=ax1) # add subplot for m value normal distribution
    ax3=fig.add_subplot(gs[1, 0], sharex=ax1) # add subplot for interactions in the range of the m plot
    ax0.set_xlim(0, chrom.cycle)
    ax1.set_xlim([0, chrom.first_95_m*2])
    
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
    global A_alpha
    A_alpha   = max(1/chrom.genesA, 0.05)
    global B_alpha
    B_alpha   = max(1/chrom.genesB, 0.05)
    global all_alpha
    all_alpha = max(1/(chrom.genesA + chrom.genesB), 0.05)
    
    # set plot titles and axis labels
    plot_title=r"$|A|={Asize}, |B|={Bsize}$, {cycles} inversion cycles".format(Asize=chrom.genesA, Bsize=chrom.genesB, cycles=chrom.cycle)
    subplot0_title=r"number of unique interactions after $n$ inversion cycles"
    subplot1_title=r"$m$ after $n$ inversion cycles"
    subplot3_title=r"number of unique interactions after $n$ inversion cycles (until cycle {})".format(chrom.first_95_m*2)
    
    ax0.set_xlabel("inversion cycle", fontsize=text_size)
    ax0.set_ylabel("unique interactions", fontsize=text_size)
    ax1.set_xlabel("inversion cycle", fontsize=text_size)
    ax1.set_ylabel(r"$m$", fontsize=text_size)
    fig.suptitle(plot_title, fontsize=plot_title_size)
    ax0.set_title(subplot0_title, fontsize=subplot_title_size)
    ax1.set_title(subplot1_title, fontsize=subplot_title_size)
    #ax3.set_title(subplot3_title, fontsize=subplot_title_size)

    return fig, ax0, ax1, ax2, ax3
    
def plot_trace(chrom, trace, ax, lim, color, alpha):
    """
    plot trace on ax from 0 to lim
    """
    for k in trace:
        ax.plot([x*chrom.sample_rate for x in range(lim)],
                 trace[k][:lim], color=color, lw = setlw, alpha=alpha)
    ax.plot([x*chrom.sample_rate for x in range(lim)],
            chrom._median_of_trace(trace)[:lim], color=color, lw = setlw*10, alpha=0.75)

def plot_t50(chrom, ax):
    """
    plot a horizontal line at x=t50 on ax
    """
    t50_text=r"$\tau_{{50\%}}={t50}$" "\n" "$({perc:.2f}\%\ of\ cycles)$".format(t50=chrom.t50, perc=chrom.t50/chrom.cycle*100)
    ax.text(x=chrom.t50, y=max(chrom._median_of_trace(chrom.trace))//10, ha='left', va='center', s=t50_text, bbox=bbox, fontsize=text_size)
    ax.axvline(x=chrom.t50, lw=setlw*5, color='black')

def plot_m(chrom, ax_m, ax_norm):
    """
    plot the m value curve and the normal distribution (on a separate plot ax_norm)
    """
    cycles=[x for x in chrom.trace_m.keys()]
    m_values=[y for y in chrom.trace_m.values()]
    pdf_space=np.linspace(0, max(m_values), 100)
    normpdf=stats.norm.pdf(pdf_space, chrom.m_mu, chrom.m_sigma)
    normpdf/=max(normpdf)
    upper_bound=chrom.m_mu+1.96*chrom.m_sigma
    lower_bound=chrom.m_mu-1.96*chrom.m_sigma
        
    ax_m.plot(cycles[:chrom.first_95_m*2+1], m_values[:chrom.first_95_m*2+1], lw=setlw*2, color='blue', label=r"$m$")
    
    # plot 95 percentile of m value normal distribution
    crossed_text="first 95 percentile value:\n{cross} cycles\n({perc:.2f}% of t50)".format(cross=chrom.first_95_m, perc=chrom.first_95_m/chrom.t50*100)
    norm_label=r"normal distribution of $m$" "\n" "(excluding the first {perc}% of cycles)".format(perc=25)
    ax_m.axhline(y=upper_bound, color='red', lw=setlw*5, ls=':') # plot upper bound of the 95 percentile
    ax_m.axhline(y=lower_bound, color='red', lw=setlw*5, ls=':') # plot lower bound of the 95 percentile
    ax_m.text(y=upper_bound, x=chrom.first_95_m*2, ha='right', va='bottom', s=r"$\mu+1.96\cdot\sigma$", bbox=bbox, fontsize=text_size)
    ax_m.text(y=lower_bound, x=chrom.first_95_m*2, ha='right', va='top', s=r"$\mu-1.96\cdot\sigma$", bbox=bbox, fontsize=text_size)
    ax_m.fill_between(cycles, lower_bound, upper_bound, color='red', alpha=0.1) # shade area between the bounds of the 95 percentile
    ax_m.axvline(x=chrom.first_95_m, lw=setlw*5, color='black') # plot the x value where the m value first enters the 95 percentile
    ax_m.text(x=chrom.first_95_m, y=0.4, ha='left', va='center', s=crossed_text, bbox=bbox, fontsize=text_size)
    
    # plot normal distribution next to m plot
    ax_norm.plot(normpdf, pdf_space, lw=setlw*5, color='red', label=norm_label) # plot the normal distribution of the m values along the y axis

def save_fig(output_dir, output_name, yaml=False):
    """
    save the figure to a specified output location as .png and .pdf
    if yaml is True, a yaml file is dumped as well
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
    
    if not yaml:
        return
    
    # save the trace as a yaml file
    if not os.path.exists(diagram_dir+'/yaml'):
        os.makedirs(diagram_dir+'/yaml')
        with open(diagram_dir+'yaml/'+output_name+'.yaml', 'w') as f:
            yaml.dump(chrom.trace, f)

def plot_results(chrom, output_dir, yaml=False):
    """
    plot the results of a simulated chromosome
    """
    # Use matplotlib to plot each key's list as a line.
    # The index of the list is the x-axis, the value is the y-axis.

    output_name=get_output_name(chrom)
    
    m_lim=(chrom.first_95_m*2+chrom.sample_rate)//chrom.sample_rate
    lim=(chrom.cycle+chrom.sample_rate)//chrom.sample_rate
    
    fig, ax0, ax1, ax2, ax3=set_up_fig(chrom)

    # plot traces
    plot_trace(chrom, chrom.trace, ax0, lim, 'black', all_alpha)
    plot_trace(chrom, chrom.trace, ax3, m_lim, 'black', all_alpha)
    plot_trace(chrom, chrom.trace_AtoB, ax0, lim, 'red', A_alpha)
    plot_trace(chrom, chrom.trace_AtoB, ax3, m_lim, 'red', A_alpha)
    plot_trace(chrom, chrom.trace_BtoA, ax0, lim, 'blue', B_alpha)
    plot_trace(chrom, chrom.trace_BtoA, ax3, m_lim, 'blue', B_alpha)  

    # plot the rest
    plot_t50(chrom, ax0)
    plot_m(chrom, ax1, ax2)

    # set legends
    ax1.legend(facecolor=fc, framealpha=box_alpha, edgecolor=ec)#, fontsize=text_size)
    ax2.legend(facecolor=fc, framealpha=box_alpha, edgecolor=ec)#, fontsize=text_size)

    save_fig(output_dir, output_name, yaml)

def log(chrom, output_dir, elapsed='N/A'):
    """
    log the state of a chromosome
    """

    # write the state to a log file according to the input
    log_dir=output_dir+'log/'
    log_file=get_output_name(chrom)+'.csv'

    if not os.path.exists(log_dir):
        os.makedirs(log_dir)

    new_log_file=not os.path.exists(log_dir+log_file)
    mode='w' if new_log_file else 'a'
    
    log_header='timestamp;|A|;|B|;cycles;t50;AB_convergence;first_95_m;m_sigma;m_mu;Delta_t\n'
    log_format_string='{ts};{A};{B};{c};{t50};{ABconv};{m95:.3f};{sig:.3f};{mu:.3f};{dt}\n'
    log_line=log_format_string.format(ts=str(dt.now()), A=chrom.genesA, B=chrom.genesB, c= chrom.cycle, t50=chrom.t50, ABconv=chrom.AB_convergence, m95=chrom.first_95_m, sig=chrom.m_sigma, mu=chrom.m_mu, dt=elapsed)

    with open(log_dir+log_file, mode) as f:
        if new_log_file:
            f.write(log_header)
        f.write(log_line)

    # write some information about the input paramters into a separate log file
    metalog_file='metalog.csv'        
    new_metalog_file=not os.path.exists(log_dir+metalog_file)

    metalog_header='|A|;|B|;converging;window_size;level_of_convergence;runs;average_t50\n'
    metalog_format_string='{A};{B};{conv};{wsize};{loc};'
    metainfo_format_string='{runs};{avt50:.2f}'
    metalog_line=metalog_format_string.format(A=chrom.genesA, B=chrom.genesB, conv=chrom.until_converged, wsize=chrom.window_size, loc=chrom.level_of_convergence)

    lines=[]
    
    if not new_metalog_file:
        with open(log_dir+metalog_file, 'r') as f:
            lines=[line.rstrip() for line in f.readlines()]
    
    with open(log_dir+metalog_file, 'w') as f:
        if new_metalog_file:
            f.write(metalog_header)
            metainfo_string=metainfo_format_string.format(runs=1, avt50=chrom.t50)
            f.write(metalog_line+metainfo_string)
        else:
            exists=False
            for line in lines:
                l=line.split(';')
                if line.startswith(metalog_line):
                    run_count=int(l[metalog_header.split(';').index('runs')])
                    run_count+=1
                    metainfo_string=metainfo_format_string.format(runs=run_count, avt50=calculate_average_t50(chrom, output_dir))
                    metalog_line+=metainfo_string
                    lines[lines.index(line)]=metalog_line
                    exists=True
                    break
            if not exists:
                metainfo_string=metainfo_format_string.format(runs=1, avt50=chrom.t50)
                lines.append(metalog_line+metainfo_string)

            for line in lines:
                f.write(line+'\n')
        
def calculate_average_t50(chrom, path):
    """
    calculate the average time (in cycles) it took to reach t50 given a set of parameters and return it
    """
    filename=get_output_name(chrom)
    raw_data=read_log_file(path, filename, cols=['t50'])
    data=[int(line[0]) for line in raw_data]
    average=np.average(data)
    return average

def plot_average_t50s(path):
    """
    read average t50 values from the metalog file and plot them
    """
    raw_data=read_log_file(path, 'metalog', cols=['|A|', '|B|', 'average_t50'])
    data=[[int(line[0])+int(line[1]), float(line[2])] for line in raw_data]

    fig, ax=plt.subplots()
    ax.scatter([line[0] for line in data], np.log10([line[1] for line in data]))
    for line in data:
        ax.annotate(r'$'+str(line[1])+'$', (line[0], np.log10(line[1])))
    ax.set_xlabel(r'$|A|+|B|$')
    ax.set_ylabel(r'$log(\bar{\tau_{50\%}})$')
    plt.savefig(path+'log/average_t50s.png')

def get_output_name(chrom):
    """
    return a string that represents the name for any output files depending on chromosome parameters (barring file endings)
    """
    return 'inversion_sim_A{Asize}-B{Bsize}_l{loc}_w{wsize}'.format(Asize=chrom.genesA, Bsize=chrom.genesB, loc=chrom.level_of_convergence, wsize=chrom.window_size)+('_c' if chrom.until_converged else '')
    
def create_parser():
    """
    creates an argparse parser for the CLI arguments
    """
    parser=ap.ArgumentParser(prog="inversion_sim", description="This program simulates inversion events of a chromosome made up of A and B genes")
    parser.add_argument('Asize', type=int, default=-1, help="integer value for the number of genes in group A")
    parser.add_argument('Bsize', type=int, default=-1, help="integer value for the number of genes in group B")
    parser.add_argument('-o', '--output-dir', default='./', help="directory in which to store the output of the program (default: './')")
    parser.add_argument('-c', '--cycle-number', type=int, default=-1, help="integer value for the number of cycles to run (alternative to --converge)")
    parser.add_argument('-C', '--converge', action='store_true', default=False, help="specify whether the simulation should run until convergence (default: False)")
    parser.add_argument('-l', '--level-of-convergence', type=float, metavar='LOC', choices=FloatRange(0, 10), default=1, help="fraction of possible gene interactions to wait for if converging (default: 1)")
    parser.add_argument('-w', '--window-size', type=int, default=1, help="the size of the window to the left and right of each gene to count as interaction after each cycle (default: 1)")
    parser.add_argument('-T', '--plot-average-t50', action='store_true', help="plot the average t50 of all previous simulations in this directory, then exit")
    
    return parser
