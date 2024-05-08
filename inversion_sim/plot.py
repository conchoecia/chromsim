import matplotlib.pyplot as plt
import matplotlib.animation as ani
import matplotlib.scale as scl
import matplotlib.patches as ptch
import matplotlib as mpl
from scipy import stats
from scipy import optimize as opt
import numpy as np
import pandas as pd
import os
from pathlib import Path

from chromosome import Chrom
import utils

def init_plot_style_settings(chrom):
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
    global marker_size
    marker_size=1000/chrom.size
    global A_alpha
    A_alpha   = max(1/chrom.Asize, 0.1)
    global B_alpha
    B_alpha   = max(1/chrom.Bsize, 0.1)
    global all_alpha
    all_alpha = max(1/(chrom.Asize + chrom.Bsize), 0.05)
    
    
def set_up_dot_fig(chrom):
    """
    set up all the axes needed for the dotplots in a grid
    """
    
    figsize=(20, 20)
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

def plot_minv(outdir, outname):
    """
    plot a boxplot of the average tS and t50, respectively
    """

    sections=utils.parse_minv_file(outdir+outname)
    
    tS_stats={k[1]: [] for k in sections.keys()}
    t50_stats={k[1]: [] for k in sections.keys()}
    for section_key in sorted(sections.keys()):
        section=sections[section_key]
        tS_stats[section_key[1]].append({'med': section['tS_median'],
                                  'q1': section['tS_q1'],
                                  'q3': section['tS_q3'],
                                  'whislo': section['tS_min'],
                                  'whishi': section['tS_max'],
                                  'label': str(section_key[0])+" genes"})
        t50_stats[section_key[1]].append({'med': section['t50_median'],
                                   'q1': section['t50_q1'],
                                   'q3': section['t50_q3'],
                                   'whislo': section['t50_min'],
                                   'whishi': section['t50_max'],
                                   'label': str(section_key[0])+" genes"})

    figsize=(20, 10)
    fig=plt.figure(figsize=figsize)
    section_count=len(sections)
    windows=sorted(list(set([section_key[1] for section_key in sections.keys()])))
    chromsizes=[sorted([section_key[0] for section_key in sections.keys() if section_key[1] == window]) for window in windows]
    
    gs=fig.add_gridspec(1, 2)
    
    ax0=fig.add_subplot(gs[0, 0])
    ax0.set_ylabel(r"$\tau_S$")
    ax0.set_xlabel("chromosome size")
    ax0.set_yscale('log')
            
    ax1=fig.add_subplot(gs[0, 1])
    ax1.set_ylabel(r"$\tau_{50}$")
    ax1.set_xlabel("chromosome size")
    ax1.set_yscale('log')

    handlers=[]
    labels=[]
    
    for i in range(0, len(windows)):
        yvals=[stats['med'] for stats in tS_stats[windows[i]]]
        color='C'+str(i)
        handlers.append(ptch.FancyBboxPatch((1, 1), 1, 1, color=color))
        labels.append("ws="+str(windows[i]))
        boxprops=dict(color=color)
        medianprops=dict(color=color)
        ax0.plot(chromsizes[i], yvals, 'D', color=color, markeredgecolor='black', linestyle='-')
        for tick in ax0.get_xticklabels():
            tick.set_rotation(-45)
        
        yvals=[stats['med'] for stats in t50_stats[windows[i]]]
        ax1.plot(chromsizes[i], yvals, 'D', color=color, markeredgecolor='black', linestyle='-')
        for tick in ax1.get_xticklabels():
            tick.set_rotation(-45)

    ax0.legend(handlers, labels)
    ax1.legend(handlers, labels)
    fig.suptitle(r"average $\tau_S$ and $\tau_{50}$")

    save_fig(outdir, outname)

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

    ax.scatter(xA_indices, yA_indices, color='blue', lw=setlw*2, s=marker_size)
    ax.scatter(xA_indices, [0]*len(xA_indices), color='blue', lw=setlw*2, s=marker_size*.5, marker='s')
    ax.scatter([0]*len(yA_indices), yA_indices, color='blue', lw=setlw*2, s=marker_size*.5, marker='s')
    ax.scatter(xB_indices, yB_indices, color='red', lw=setlw*2, s=marker_size)
    ax.scatter(xB_indices, [0]*len(xB_indices), color='red', lw=setlw*2, s=marker_size*.5, marker='s')
    ax.scatter([0]*len(yB_indices), yB_indices, color='red', lw=setlw*2, s=marker_size*.5, marker='s')
    
def plot_trace(chrom, trace, lim, ax, mean, color='blue', alpha=0.5):
    """
    plot trace dictionary on the given axes in the interval [0, lim]
    """

    # plot all the individual traces
    for k in trace:
        x, y= [], []
        for cycle in trace[k]:
            if cycle >= lim:
                break
            x.append(cycle)
            y.append(trace[k][cycle])
        ax.plot(x, y, color=color, lw = setlw, alpha=alpha)

    # get the x and y lists from mean
    x, y=[], []
    for cycle in mean:
        if cycle >= lim:
            break
        x.append(cycle)
        y.append(mean[cycle])
        
    # plot the fitted data
    ax.plot(x, y, color=color, lw = setlw*10, alpha=0.75)

def plot_t50(chrom, ax):
    """
    plot a horizontal line at x=t50 on the given axis
    """
    
    t50_text=r"$\tau_{{50\%}}={t50}$" "\n" r"$({perc:.2f}\%\ of\ cycles)$".format(t50=chrom.t50, perc=chrom.t50/chrom.cycle*100)
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
                
    # save this as a pdf and png
    plt.savefig(outdir+outname+'.pdf')
    plt.savefig(outdir+outname+'.png')
    
def plot_chrom(source, outdir, gif=False):
    """
    plot a chromosome from a .inv source
    """
    
    chrom, results, ts=utils.parse_inv_file(source)

    # get the file name without the file ending
    outname=Path(source).stem
    
    init_plot_style_settings(chrom)
    
    if gif:
        print("generating GIF file of the inversion process - this might take a while")
        
        make_dotplot_gif(chrom, outdir, results['tS'], outname)
    else:
        chrom.run(len(chrom.inversion_cuts))

        plot_results(chrom, outdir, outname)

        
def make_dotplot_gif(chrom, outdir, cycles, outname):
    """
    create an animated GIF of the chromosome's inversion process up to cycles (max. 1000 frames)
    """

    max_frames=1000
    frames=min(max_frames, cycles)
    frames=cycles
    step=cycles/frames

    figsize=(10, 10) if chrom.size <= 400 else (5, 5)
    fig=plt.figure(figsize=figsize)
    global marker_size
    marker_size/=8
    ax=fig.add_subplot()

    last_seen=0
    def animate(i):
        chrom.run(step*i)
        ax.clear()
        nonlocal last_seen
        ax.text(0, 1.05, r"cycle: {c:4d}, m: {m:.2f}, seen (new): {n}, seen (cumulative): {a}".format(c=chrom.cycle, m=chrom.trace_m[chrom.cycle], n=len(chrom.seen)-last_seen, a=len(chrom.seen)), transform=ax.transAxes, ha='left', weight='bold')
        last_seen=len(chrom.seen)
        plot_dotplot(chrom, ax, chrom.original_gene_list, chrom.gene_list)
        
    animation=ani.FuncAnimation(fig, animate, repeat=True, frames=frames, interval=500)
    
    animation.save(outdir+outname+'.gif', writer='imagemagick')

def plot_results(chrom, outdir, outname):
    """
    plot the results of a simulated chromosome
    """
    
    ## trace figure

    # calculate the upper bounds up to which should be plotted
    m_lim=chrom.tS*2+1
    lim=len(chrom.inversion_cuts)

    # set up the fig
    fig, ax0, ax1, ax2, ax3=set_up_trace_fig(chrom)

    # plot traces
    plot_trace(chrom, chrom.trace, lim, ax0, chrom.trace_mean, 'black', all_alpha)
    plot_trace(chrom, chrom.trace, m_lim, ax3, chrom.trace_mean, 'black', all_alpha)
    plot_trace(chrom, chrom.trace_AtoB, lim, ax0, chrom.trace_AtoB_mean, 'red', A_alpha)
    plot_trace(chrom, chrom.trace_AtoB, m_lim, ax3, chrom.trace_AtoB_mean, 'red', A_alpha)
    plot_trace(chrom, chrom.trace_BtoA, lim, ax0, chrom.trace_BtoA_mean, 'blue', B_alpha)
    plot_trace(chrom, chrom.trace_BtoA, m_lim, ax3, chrom.trace_BtoA_mean, 'blue', B_alpha)  

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
    ts_chrom=Chrom(chrom.Asize, chrom.Bsize, level_of_convergence=chrom.level_of_convergence, window_size=chrom.window_size, inversion_cuts=chrom.inversion_cuts)
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
