import matplotlib.pyplot as plt
from scipy import stats
import numpy as np
import pandas as pd
import os

from chromosome import Chrom
import utils

def get_output_name(chrom):
    """
    return a string to be used as the base for saving logs

    TODO: will be redundant with the new log file format
    """
    
    return 'inversion_sim_A{Asize}-B{Bsize}_l{loc}_w{wsize}'.format(Asize=chrom.genesA, Bsize=chrom.Bsize, loc=chrom.level_of_convergence, wsize=chrom.window_size)

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

def set_up_dot_fig(chrom):
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

    ax0.set_title(subplot0_title, fontsize=subplot_title_size)
    ax1.set_title(subplot1_title, fontsize=subplot_title_size)
    ax2.set_title(subplot2_title, fontsize=subplot_title_size)
    ax3.set_title(subplot3_title, fontsize=subplot_title_size)

    return fig, ax0, ax1, ax2, ax3
    
def set_up_trace_fig(chrom):
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
    A_alpha   = max(1/chrom.Asize, 0.1)
    global B_alpha
    B_alpha   = max(1/chrom.Bsize, 0.1)
    global all_alpha
    all_alpha = max(1/(chrom.Asize + chrom.Bsize), 0.05)
    
    # set plot titles and axis labels
    subplot0_title=r"number of unique interactions after $n$ inversion cycles"
    subplot1_title=r"$m$ after $n$ inversion cycles"

    ax0.set_xlabel("inversion cycle", fontsize=text_size)
    ax0.set_ylabel("unique interactions", fontsize=text_size)
    ax1.set_xlabel("inversion cycle", fontsize=text_size)
    ax1.set_ylabel(r"$m$", fontsize=text_size)
    ax3.set_xlabel("inversion cycle", fontsize=text_size)
    ax3.set_ylabel("unique interactions", fontsize=text_size)

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

def save_fig(outdir, outname):
    """
    save the figure of the current matplotlib plot to a specified output location as .png and .pdf
    """
    
    # create the diagram directory if it does not exist yet
    diagram_dir=outdir+'diagrams/'
    
    if not os.path.exists(diagram_dir):
        os.makedirs(diagram_dir)
        if not os.path.exists(diagram_dir+'/pdf'):
            os.makedirs(diagram_dir+'/pdf')
            if not os.path.exists(diagram_dir+'/png'):
                os.makedirs(diagram_dir+'/png')
                
    # save this as a pdf and png
    plt.savefig(outname+'.pdf')
    plt.savefig(outname+'.png')

def plot_chrom(source, outdir, gif=False):
    chrom, results, ts=utils.parse_inv_file(source)

    if not gif:
        chrom.run(len(chrom.inversion_cuts))

        plot_results(chrom, outdir)
    
def plot_results(chrom, outdir):
    """
    plot the results of a simulated chromosome
    """    

    init_plot_style_settings()

    outname=str(chrom.timestamp)
    
    ## trace figure

    # calculate the upper bounds up to which should be plotted
    m_lim=(chrom.tS*2)
    lim=(chrom.cycle)

    # set up the fig
    fig, ax0, ax1, ax2, ax3=set_up_trace_fig(chrom)

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
    save_fig(outdir, outname+"_trace")

    ## dotplot figure

    # get the gene lists
    original_genes=chrom.original_gene_list
    final_genes=chrom.gene_list
    half_time_genes=chrom.t50_gene_list

    # rerun the chromosome until tS to get the list at that point
    ts_chrom=Chrom(chrom.Asize, chrom.Bsize, chrom.level_of_convergence, chrom.window_size, chrom.translocations_per_cycle, chrom.inversion_cuts)
    ts_chrom.run(n=chrom.tS, show_output=False)
    entropy_genes=ts_chrom.gene_list

    # set up the fig
    fig, ax0, ax1, ax2, ax3=set_up_dot_fig(chrom)

    # plot dotplots
    plot_dotplot(chrom, ax0, original_genes, final_genes)
    if chrom.t50 >= 0:
        plot_dotplot(chrom, ax1, original_genes, half_time_genes)
    if chrom.tS >= 0:
        plot_dotplot(chrom, ax2, original_genes, entropy_genes)
    if chrom.tS >= 0 and chrom.t50 >= 0:
        plot_dotplot(chrom, ax3, entropy_genes, half_time_genes)

    # save the dotplot figure
    save_fig(outdir, outname+"_dotplots")
