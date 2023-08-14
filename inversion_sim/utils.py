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
    output_file='runtime.png'
    raw_data=read_log_file(path, filename, ['|A|', '|B|', 'average_runtime'])
    data=[[int(line[0])+int(line[1]), line[2]] for line in raw_data]
    
    plt.plot([x[0] for x in data], [runtime_to_int(x[1])/60 for x in data], 'o')
    plt.xlabel("chromosome size (|A|+|B|)")
    plt.ylabel("runtime (minutes)")
    plt.savefig(path+output_file)
    print('success')
    
def init_plot_style_settings():
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

def set_up_dot_fig(chrom):
    """
    TODO
    """
    figsize=(35, 10)
    fig=plt.figure(figsize=figsize)
    gs=fig.add_gridspec(1, 3)
    ax0=fig.add_subplot(gs[0, 0])
    ax1=fig.add_subplot(gs[0, 1])
    ax2=fig.add_subplot(gs[0, 2])

    plot_title=r"$|A|={Asize}, |B|={Bsize}$, {cycles} inversion cycles".format(Asize=chrom.genesA, Bsize=chrom.genesB, cycles=chrom.cycle)
    subplot0_title=r"dotplot start vs. $\tau{100\%}$ of simulation"
    subplot1_title=r"dotplot start vs. $\tau_{50\%}$ of simulation"
    subplot2_title=r"dotplot start vs. $\tau_{S}$ of simulation"

    ax0.set_xlabel("start", fontsize=text_size)
    ax0.set_ylabel(r"$\tau{100\%}$", fontsize=text_size)
    ax1.set_xlabel("start", fontsize=text_size)
    ax1.set_ylabel(r"$\tau{50\%}$", fontsize=text_size)
    ax2.set_xlabel("start", fontsize=text_size)
    ax2.set_ylabel(r"$\tau{S}$", fontsize=text_size)

    ax0.set_title(subplot0_title, fontsize=subplot_title_size)
    ax1.set_title(subplot1_title, fontsize=subplot_title_size)
    ax2.set_title(subplot2_title, fontsize=subplot_title_size)

    return fig, ax0, ax1, ax2
    
def set_up_trace_fig(chrom):
    """
    create the basic structure of the matplotlib figure
    returns the figure itself and the axes that are part of it
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
    plot_title=r"$|A|={Asize}, |B|={Bsize}$, {cycles} inversion cycles".format(Asize=chrom.genesA, Bsize=chrom.genesB, cycles=chrom.cycle)
    subplot0_title=r"number of unique interactions after $n$ inversion cycles"
    subplot1_title=r"$m$ after $n$ inversion cycles"
    subplot3_title=r"number of unique interactions after $n$ inversion cycles (until cycle {})".format(chrom.tS*2)
    
    ax0.set_xlabel("inversion cycle", fontsize=text_size)
    ax0.set_ylabel("unique interactions", fontsize=text_size)
    ax1.set_xlabel("inversion cycle", fontsize=text_size)
    ax1.set_ylabel(r"$m$", fontsize=text_size)
    ax3.set_xlabel("inversion cycle", fontsize=text_size)
    ax3.set_ylabel("unique interactions", fontsize=text_size)
    fig.suptitle(plot_title, fontsize=plot_title_size)
    ax0.set_title(subplot0_title, fontsize=subplot_title_size)
    ax1.set_title(subplot1_title, fontsize=subplot_title_size)
    #ax3.set_title(subplot3_title, fontsize=subplot_title_size)

    return fig, ax0, ax1, ax2, ax3

def plot_dotplot(chrom, ax, x, y):
    """
    plot a dotplot between the original gene list and the final ones
    """
    x_indices=[x.index(gene) for gene in x]
    y_indices=[y.index(gene) for gene in x]
    ax.scatter(x_indices, y_indices, lw=setlw*2, color='blue')

def plot_trace(chrom, trace, ax, lim, color, alpha):
    """
    plot trace on ax from 0 to lim
    """
    all=[]
    for k in trace:
        x, y= [], []
        for tup in trace[k]:
            if tup[0] > lim:
                break
            x.append(tup[0])
            y.append(tup[1])
            all.append((tup[0], tup[1]))
        ax.plot(x, y, color=color, lw = setlw, alpha=alpha)
    all=sorted(all)

    # calculate the moving average of the combined list
    df = pd.DataFrame(all, columns =['cycle', 'interactions'])
    new_all=df.rolling(len(all)//20).mean()

    # add the last elements of the original list to make sure the line goes all the way to the end
    x=new_all['cycle'].tolist()
    x.append(all[-1][0])
    y=new_all['interactions'].tolist()
    y.append(all[-1][1])
    
    ax.plot(x, y, color=color, lw = setlw*10, alpha=0.75)

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

def plot_results(chrom, output_dir, do_yaml=False):
    """
    plot the results of a simulated chromosome
    """    
    init_plot_style_settings()

    output_name=get_output_name(chrom)
    
    m_lim=(chrom.tS*2)#+chrom.sample_rate)//chrom.sample_rate
    lim=(chrom.cycle)#+chrom.sample_rate)//chrom.sample_rate
    
    fig, ax0, ax1, ax2, ax3=set_up_trace_fig(chrom)

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

    save_fig(output_dir, output_name+"_trace", do_yaml)

    x=chrom.original_gene_list
    y0=chrom.gene_list
    y1=chrom.t50_gene_list
    y2=chrom.tS_gene_list
    
    fig, ax0, ax1, ax2=set_up_dot_fig(chrom)

    plot_dotplot(chrom, ax0, x, y0)
    plot_dotplot(chrom, ax1, x, y1)
    #plot_dotplot(chrom, ax2, x, y2)

    save_fig(output_dir, output_name+"_dotplots")

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
    log_format_string='{ts};{A};{B};{c};{t50};{ABconv};{m95};{sig:.3f};{mu:.3f};{dt}\n'
    log_line=log_format_string.format(ts=str(dt.now()), A=chrom.genesA, B=chrom.genesB, c= chrom.cycle, t50=chrom.t50, ABconv=chrom.AB_convergence, m95=chrom.tS, sig=chrom.m_sigma, mu=chrom.m_mu, dt=elapsed)

    with open(log_dir+log_file, mode) as f:
        if new_log_file:
            f.write(log_header)
        f.write(log_line)

def runtime_to_int(runtime):
    arr=runtime.split(':')
    mins=int(arr[0])
    secs=int(arr[1])
    return mins*60+secs

def int_to_runtime(f):
    return str(int(f//60))+':'+str(int(f%60)).zfill(2)
        
def metalog(output_dir):
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
                runs=[run.rstrip().split(';') for run in lf.readlines()]

            header=runs[0]
            t50s=[float(run[header.index('t50')]) for run in runs[1:]]
            t100s=[float(run[header.index('cycles')]) for run in runs[1:]]
            tss=[float(run[header.index('first_95_m')]) for run in runs[1:]]
            rts=[runtime_to_int(run[header.index('Delta_t')]) for run in runs[1:]]
            
            metalog_line=metalog_format_string.format(A=Asize, B=Bsize, conv=converging, wsize=wsize, loc=loc)        
            metainfo_string=metainfo_format_string.format(runs=len(runs), avt50=np.average(t50s), avt100=np.average(t100s), avts=np.average(tss), avrt=int_to_runtime(np.average(rts)))
            mf.write(metalog_line+metainfo_string)
        
def calculate_average(chrom, path, col):
    """
    calculate the average time (in cycles) it took to reach t50 given a set of parameters and return it
    """
    filename=get_output_name(chrom)
    raw_data=read_log_file(path, filename, cols=[col])
    data=[int(line[0]) for line in raw_data]
    average=np.average(data)
    return average
        
#def calculate_average_t100(chrom, path):
#,average_tentrop;{avtent:.2f}y    """
#,average_tentropy    calculate the average time (in cycles) it took to reach t50 given a set of parameters and return it
#,average_tentropy    """
#,average_tentrop;{avtent:.2f}y    filename=get_output_name(chrom)
#,average_tentrop;{avtent:.2f}y    raw_data=read_log_file(path, filename, cols=['cycles'])
#,average_tentrop;{avtent:.2f}y    data=[int(line[0]) for line in raw_data]
#,average_tentrop;{avtent:.2f}y    average=np.average(data)
#,average_tentrop;{avtent:.2f}y    return average;{avtent:.2f}

def plot_average_t50s(path):
    """
    read average t50 values from the metalog file and plot them
    """
    cols=['|A|', '|B|', 'average_t50', 'average_t100', 'average_ts', 'window_size']
    raw_data=read_log_file(path, 'metalog', cols=cols)
    data=[[int(line[cols.index('|A|')])+int(line[cols.index('|B|')]),
           float(line[cols.index('average_t50')]),
           float(line[cols.index('average_t100')]),
           float(line[cols.index('average_ts')]),
           int(line[cols.index('window_size')])]
          for line in raw_data if float(line[2]) >= 0]

    wsizes=[]
    grouped_data={}
    for line in data:
        wsize=line[-1]
        if wsize not in grouped_data:
            grouped_data[wsize]=[]
        grouped_data[wsize].append([line[0], line[1], line[2], line[3]])
    
    init_plot_style_settings()
    fig=plt.figure(figsize=(20, 20))
    gs=fig.add_gridspec(4, 2)
    ax0=fig.add_subplot(gs[0, 0])
    ax1=fig.add_subplot(gs[1, 0])
    ax2=fig.add_subplot(gs[0, 1])
    ax3=fig.add_subplot(gs[1, 1])
    ax4=fig.add_subplot(gs[2, 0])
    ax5=fig.add_subplot(gs[3, 0])
    ax6=fig.add_subplot(gs[2, 1])
    ax7=fig.add_subplot(gs[3, 1])

    for k in grouped_data:
        kdata=grouped_data[k]
        x=[line[0] for line in kdata]
        y=[line[1] for line in kdata]
        logy=np.log(y)
        ydiff=[line[2]-line[1] for line in kdata]
        logydiff=np.log(ydiff)
        yratio=[line[1]/line[2] for line in kdata]
        logyratio=np.log(yratio)
        ydiffs=[line[1]-line[3] for line in kdata]
        logydiffs=np.log(ydiffs)

        label='w{}'.format(k)

        ax0.scatter(x, y, label=label)
        ax1.scatter(x, logy, label=label)
        ax2.scatter(x, ydiff, label=label)
        ax3.scatter(x, logydiff, label=label)
        ax4.scatter(x, yratio, label=label)
        ax5.scatter(x, logyratio, label=label)
        ax6.scatter(x, ydiffs, label=label)
        ax7.scatter(x, logydiffs, label=label)
    
    #for line in data:
    #    offset=max([np.log10(l[1]) for l in data])
    #    print(offset)
    #    print(offset*0.2)
    #    print(np.log10(line[1])+offset*0.2)
    #    ax0.annotate('w'+str(line[2]), (line[0], line[1]+offset*0.02), bbox=bbox, fontsize=text_size/2)
    #    ax1.annotate('w'+str(line[2]), (line[0], np.log10(line[1])+offset*0.02), bbox=bbox, fontsize=text_size/2)
    #ax0.set_xlabel(r'$|A|+|B|$')
    ax5.set_xlabel(r'$|A|+|B|$', fontsize=text_size)
    ax7.set_xlabel(r'$|A|+|B|$', fontsize=text_size)
    ax0.set_ylabel(r'$\bar{\tau}_{50\%}$', fontsize=text_size)
    ax1.set_ylabel(r'$ln(\bar{\tau}_{50\%})$', fontsize=text_size)
    ax2.set_ylabel(r'$\bar{\tau}_{50\%}-\bar{\tau}_{100\%}$', fontsize=text_size)
    ax3.set_ylabel(r'$ln(\bar{\tau}_{50\%}-\bar{\tau}_{100\%})$', fontsize=text_size)
    ax4.set_ylabel(r'$\frac{\bar{\tau}_{50\%}}{\bar{\tau}_{100\%}}$', fontsize=text_size)
    ax5.set_ylabel(r'$ln(\frac{\bar{\tau}_{50\%}}{\bar{\tau}_{100\%}})$', fontsize=text_size)
    ax6.set_ylabel(r'$\bar{\tau}_{50\%}-\bar{\tau}_{S}$', fontsize=text_size)
    ax7.set_ylabel(r'$ln(\bar{\tau}_{50\%}-\bar{\tau}_{S})$', fontsize=text_size)
    ax0.legend()
    ax1.legend()
    ax2.legend()
    ax3.legend()
    ax4.legend()
    ax5.legend()
    ax6.legend()
    ax7.legend()

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
    # flag simulation  arguments
    parser.add_argument('-C', '--converge', action='store_true', default=False, help="specify whether the simulation should run until convergence (default: False)")
    parser.add_argument('-P', '--plot-curves', action='store_true', default=False, help="tell the simulation to plot the curves in the end (default: False)")
    # flag simulation-free arguments
    parser.add_argument('-M', '--metalog', action='store_true', help="collect averages from log files and store them in metalog.csv (not running a simulation)")
    parser.add_argument('-T', '--plot-average-t50', action='store_true', help="plot the average t50 of all previous simulations in this directory (not running a simulation)")
    
    return parser
