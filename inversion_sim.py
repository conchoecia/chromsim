#!/usr/bin/env python

"""
This program calculates the time it takes for any gene to interact with another
after a fusion-with mixing event.
"""

import random
import sys
import time
import os

remote_path_base='/scratch/molevo/bluehmel/'

class Chrom():
    def __init__(self, size, gene_quantityA, gene_quantityB):
        self.length = size # is this even relevant when gene_quantityA and gene_quantityB are specified?
        self.genesA = gene_quantityA
        self.genesB = gene_quantityB
        self.size=self.genesA+self.genesB

        self.gene_list = ["A."+str(i+1) for i in range(self.genesA)] + ["B."+str(i+1) for i in range(self.genesB)]
        # seen is a list of genes that have already interacted with each other, value is cycle
        self.seen = {}
        self.t50=-1
        self.AB_convergence=-1
        #self.B_convergence=-1

        """if |A| + |B| <= 50, sample rate 1.0
elif |A| + |B| <= 100, sample rate 0.5
elif |A| + |B| <= 500, sample rate 0.1
else sample rate 0.01"""
        # this is the percent of the data that we sample, range [0-1]. 0.01 is a good rate for samples above 100k iterations
        self.sample_frequency = 1 if self.size <= 50 else 0.5 if self.size <= 100 else 0.1 if self.size <= 500 else 0.01
        self.sample_rate = int( 1 / self.sample_frequency)

        # this tracks how many new interactions were added at each cycle
        self.trace      = {k:[] for k in self.gene_list}
        self.trace_AtoB = {k:[] for k in self.gene_list if k.startswith("A")}
        self.trace_BtoA = {k:[] for k in self.gene_list if k.startswith("B")}
        self.trace_m={}
        self.cycle = 0
        self.update_seen()
        self.calculate_m()
        self.first_95_m=-1
        self.m_sigma=-1
        self.m_mu=-1

        # the following is just for logging/debugging purposes
        self.last_i0=0
        self.last_i1=self.genesA+self.genesB
        #self.log=[str(self)]

    def __str__(self):
        format_string="cycle: {cycle:10d}; last inversion: {start} {istart:5d} {inverted} {iend:<5d} {end}; m: {m:1.3f}"
        AB_string=self.get_AB_string()
        i0=self.last_i0
        i1=self.last_i1
        m=self.trace_m[self.cycle]
        ret=format_string.format(cycle=self.cycle, start=AB_string[0:i0], istart=i0, inverted=AB_string[i0:i1], iend=i1-1, end=AB_string[i1:len(AB_string)], m=m)
        return ret
    
    def simulation_cycle(self, iterations = 0, until_converged = False):
        """
        This function runs the simulation for a given number of iterations, or until done.
        """
        # raise an error if we don't know how to run this method
        if iterations == 0 and until_converged == False:
            raise ValueError("Please specify either iterations or until_converged.")
        # until_converged hasn't been implemented yet, so tell the user that this isn't possible now
        if until_converged:
            converging_at=self.calculate_convergence()
            AB_have_converged=False
            #B_has_converged=False
            level_of_convergence=.99 # the percentage of possible interactions to be waiting for before stopping to limit runtime to a feasible range
            while len(self.seen)/converging_at < level_of_convergence: 
                self.shuffle()
                if self.t50 < 0 and len(self.seen) >= converging_at/2:
                    self.t50=self.cycle
                    print("reached t50 at "+str(self.t50))
                if self.t50 >= 0:
                    if self.AB_convergence < 0:
                        AB_have_converged=True
                        for k in self.trace_AtoB:
                            if self.trace_AtoB[k][-1] < self.genesB-(1 if self.gene_list.index(k) == 0 else 0):
                                AB_have_converged=False
                    #if self.B_convergence < 0:
                    #    B_has_converged=True
                    #    for k in self.trace_BtoA:
                    #        if self.trace_BtoA[k][-1] < self.genesA-(1 if self.gene_list.index(k) == self.size-1 else 0):
                    #            B_has_converged=False
                    if AB_have_converged and self.AB_convergence < 0:
                        self.AB_convergence=self.cycle
                        print("A-B/B-A converged at "+str(self.AB_convergence))
                    #if B_has_converged and self.B_convergence < 0:
                    #    self.B_convergence=self.cycle
                    #    print("B converged at "+str(self.B_convergence))
            print("converged at "+str(self.cycle))
        else:
            # run the simulation for the specified number of iterations
            for i in range(iterations):
                # make a progress bar that stays on the same line that prints every 100 cycles
                if i % 100 == 0:
                    # use a carriage return to stay on the same line
                    # use format string to make the number occupy 10 spaces.
                    # separate with a space character, then plot the percentage.
                    print("\r{cycle:15d} {progress:.2f}%  ".format(Asize=self.genesA, Bsize=self.genesB, cycle=i, progress=(i/iterations)*100), end="")
                self.shuffle()
            print("{cycle:15d}  100.00%".format(cycle=i+1))

    """
    calculates the number of possible unique gene interactions
    """
    def calculate_convergence(self):
        conv=0
        for i in range(self.genesA+self.genesB):
            conv+=i
        return conv-1
    
    def shuffle(self):
        # Randomly pick two indices in the list.
        # Start at one and end at len-1 to not destroy telomeres
        i0 = random.randint(1, len(self.gene_list)-1)
        i1 = random.randint(1, len(self.gene_list)-1)
        sortedi = sorted([i0, i1]) 
        i0 = sortedi[0]
        i1 = sortedi[1]
        self.gene_list[i0:i1] = self.gene_list[i0:i1][::-1]
        self.update_seen()
        self.cycle += 1
        self.calculate_m()

        #log the current state
        self.last_i0=i0
        self.last_i1=i1
        #self.log.append(str(self))
    
    def update_seen(self):
        """
        update the seen graph
        """
        # increment the trace structure so we can modify
        for k in self.trace:
            if len(self.trace[k]) == 0:
                self.trace[k].append(0)
                if k in self.trace_AtoB:
                    self.trace_AtoB[k].append(0)
                if k in self.trace_BtoA:
                    self.trace_BtoA[k].append(0)
            elif self.cycle % self.sample_rate == 0:
                self.trace[k].append(self.trace[k][-1])
                if k in self.trace_AtoB:
                    self.trace_AtoB[k].append(self.trace_AtoB[k][-1])
                if k in self.trace_BtoA:
                    self.trace_BtoA[k].append(self.trace_BtoA[k][-1])
        # go through all of the pairs in the list to update self.seen
        for i in range(len(self.gene_list)-1):
            this_edge = tuple(sorted([self.gene_list[i],
                                      self.gene_list[i+1]]))
            if this_edge not in self.seen:
                #if "B.1" in this_edge:
                #    print(this_edge)
                #if "A.1" in this_edge:
                #    print(this_edge)
                # add the edge to the seen graph
                self.seen[this_edge] = self.cycle
                # update the trace structure
                for j in [0,1]:
                    other = 1 if j == 0 else 0
                    self.trace[this_edge[j]][-1] += 1
                    if this_edge[j].startswith("A") and this_edge[other].startswith("B"):
                        self.trace_AtoB[this_edge[j]][-1] += 1
                    if this_edge[j].startswith("B") and this_edge[other].startswith("A"):
                        self.trace_BtoA[this_edge[j]][-1] += 1

    def get_AB_string(self):
        return ''.join([gene[0] for gene in self.gene_list])

    """
    calculate m of the current gene sequence
    """
    def calculate_m(self):
        sequence=self.get_AB_string()
        substrings = [sequence[i:i+2] for i in range(len(sequence)-1)]
        A = sequence.count('A')
        B = sequence.count('B')
        AB = substrings.count('AB')
        BA = substrings.count('BA')
        m = (AB + BA - 1)/ ((2* A * B)/(A+B) - 1)
        self.trace_m[self.cycle]=m
        
    def _median(self, lst):
        sortedLst = sorted(lst)
        lstLen = len(lst)
        index = (lstLen - 1) // 2
       
        if (lstLen % 2):
            return sortedLst[index]
        else:
            return (sortedLst[index] + sortedLst[index + 1])/2.0
    
    def _median_of_trace(self, trace):
        """
        Gets the median value of the supplied traces.
        """
        # get a random key of self.trace
        k = list(trace.keys())[0]
        # calculate the median of all the traces at each sampling point
        return [self._median([trace[j][i] for j in trace]) for i in range(len(trace[k]))]

    def plot_results(self, yaml=False):
        """
        Plot the trace as faint lines.
        """
        # Use matplotlib to plot each key's list as a line.
        # The index of the list is the x-axis, the value is the y-axis.

        # calculate all the values
        from scipy import stats
        import math
        import numpy as np

        cycles=[x for x in self.trace_m.keys()]
        m_values=[y for y in self.trace_m.values()]
        burn_in=0.25
        start_norm_at=math.floor((self.cycle+1)*burn_in)
        burnt_m_values=np.array(m_values[start_norm_at:])
        mu=np.mean(burnt_m_values)
        var=np.var(burnt_m_values)
        sigma=math.sqrt(var)
        self.m_sigma=sigma
        self.m_mu=mu
        pdf_space=np.linspace(0, max(burnt_m_values), 100)
        normpdf=stats.norm.pdf(pdf_space, mu, sigma)
        normpdf/=max(normpdf)
        upper_bound=mu+1.96*sigma
        lower_bound=mu-1.96*sigma
        crossed_lower_bound_at=0
        for k in self.trace_m.items():
            if k[1] >= lower_bound:
                crossed_lower_bound_at=k[0]
                break
        self.first_95_m=crossed_lower_bound_at
        cycle_limit=crossed_lower_bound_at*2+1
        sample_limit=cycle_limit//self.sample_rate+1

        # initialize the matplotlib plot
        import matplotlib.pyplot as plt

        #print([k for k in  self.trace_AtoB])
        #print("A.1: ", self.trace_AtoB["A.1"])
        #print([k for k in  self.trace_BtoA], self.trace_BtoA)

        # set up the figure
        figsize=(20, 25)
        fig=plt.figure(figsize=figsize)
        gs=fig.add_gridspec(3, 2, width_ratios=[3, 1])
        ax0=fig.add_subplot(gs[0, :]) # add subplot for gene interaction traces spanning the top half
        ax1=fig.add_subplot(gs[2, 0]) # add subplot for m values
        ax2=fig.add_subplot(gs[2, 1], sharey=ax1) # add subplot for m value normal distribution
        ax3=fig.add_subplot(gs[1, 0], sharex=ax1) # add subplot for interactions in the range of the m plot
        ax0.set_xlim(0, self.cycle)
        ax1.set_xlim([0, cycle_limit])

        # figure styling
        plt.style.use('bmh')
        setlw = 0.2
        seta  = 0.1
        A_alpha   = max(1/self.genesA, 0.05)
        B_alpha   = max(1/self.genesB, 0.05)
        all_alpha = max(1/(self.genesA + self.genesB), 0.05)
        fc='white'
        box_alpha=0.5
        ec='gray'
        bbox=dict(facecolor=fc, alpha=box_alpha, boxstyle='round, pad=0.5', edgecolor=ec)
        plot_title_size=50
        subplot_title_size=30
        text_size=20

        # set plot titles
        plot_title=r"$|A|={Asize}, |B|={Bsize}$, {cycles} inversion cycles".format(Asize=self.genesA, Bsize=self.genesB, cycles=self.cycle)
        subplot0_title=r"number of unique interactions after $n$ inversion cycles"
        subplot1_title=r"$m$ after $n$ inversion cycles"
        subplot3_title=r"number of unique interactions after $n$ inversion cycles (until cycle {})".format(cycle_limit)

        # plot gene interaction trace
        for k in self.trace:
            ax0.plot([x*self.sample_rate for x in range(len(self.trace[k]))],
                     self.trace[k], color='black', lw = setlw, alpha=all_alpha)
            ax3.plot([x*self.sample_rate for x in range(cycle_limit+1)],
                     self.trace[k][:cycle_limit+1], color='black', lw = setlw, alpha=all_alpha)
        median_of_trace=self._median_of_trace(self.trace)
        ax0.plot([x*self.sample_rate for x in range(len(self.trace[k]))],
                 median_of_trace, color='black', lw = setlw*10, alpha=0.75)
        ax3.plot([x*self.sample_rate for x in range(cycle_limit+1)],
                 median_of_trace[:cycle_limit+1], color='black', lw = setlw*10, alpha=0.75)
 
        # now plot the A-to-B and B-to-A traces
        for k in self.trace_AtoB:
            ax0.plot([x*self.sample_rate for x in range(len(self.trace_AtoB[k]))],
                     self.trace_AtoB[k], color='blue', lw = setlw, alpha=A_alpha)
            ax3.plot([x*self.sample_rate for x in range(cycle_limit+1)],
                     self.trace_AtoB[k][:cycle_limit+1], color='blue', lw = setlw, alpha=A_alpha)
        ax0.plot([x*self.sample_rate for x in range(len(self.trace_AtoB[k]))],
                 self._median_of_trace(self.trace_AtoB), color='blue', lw = setlw*10, alpha=0.75)
        ax3.plot([x*self.sample_rate for x in range(cycle_limit+1)],
                 self._median_of_trace(self.trace_AtoB)[:cycle_limit+1], color='blue', lw = setlw*10, alpha=0.75)

        for k in self.trace_BtoA:
            ax0.plot([x*self.sample_rate for x in range(len(self.trace_BtoA[k]))],
                     self.trace_BtoA[k], color='red', lw = setlw, alpha=B_alpha)
            ax3.plot([x*self.sample_rate for x in range(cycle_limit+1)],
                     self.trace_BtoA[k][:cycle_limit+1], color='red', lw = setlw, alpha=B_alpha)
        ax0.plot([x*self.sample_rate for x in range(len(self.trace_BtoA[k]))],
                 self._median_of_trace(self.trace_BtoA), color='red', lw = setlw*10, alpha=0.75)
        ax3.plot([x*self.sample_rate for x in range(cycle_limit+1)],
                 self._median_of_trace(self.trace_BtoA)[:cycle_limit+1], color='red', lw = setlw*10, alpha=0.75)

        t50_text=r"$\tau_{{50\%}}={t50}$" "\n" "$({perc:.2f}\%\ of\ cycles)$".format(t50=self.t50, perc=self.t50/self.cycle*100)
        ax0.text(x=self.t50, y=max(median_of_trace)//10, ha='left', va='center', s=t50_text, bbox=bbox, fontsize=text_size)
        ax0.axvline(x=self.t50, lw=setlw*5, color='black')
        
        #plot m values on the second subplot
        ax1.plot(cycles[:cycle_limit+1], m_values[:cycle_limit+1], lw=setlw*2, color='blue', label=r"$m$")

        # plot 95 percentile of m value normal distribution
        crossed_text="first 95 percentile value:\n{cross} cycles\n({perc:.2f}% of t50)".format(cross=crossed_lower_bound_at, perc=crossed_lower_bound_at/self.t50*100)
        burn_in_text="burn-in:\n{cycle} cycles\n({perc:.0f}% of cycles)".format(cycle=start_norm_at, perc=burn_in*100)
        norm_label=r"normal distribution of $m$" "\n" "(excluding the first {perc}% of cycles)".format(perc=burn_in*100)
        ax1.axhline(y=upper_bound, color='red', lw=setlw*5, ls=':') # plot upper bound of the 95 percentile
        ax1.axhline(y=lower_bound, color='red', lw=setlw*5, ls=':') # plot lower bound of the 95 percentile
        ax1.text(y=upper_bound, x=cycle_limit, ha='right', va='bottom', s=r"$\mu+1.96\cdot\sigma$", bbox=bbox, fontsize=text_size)
        ax1.text(y=lower_bound, x=cycle_limit, ha='right', va='top', s=r"$\mu-1.96\cdot\sigma$", bbox=bbox, fontsize=text_size)
        ax1.fill_between(cycles, lower_bound, upper_bound, color='red', alpha=0.1) # shade area between the bounds of the 95 percentile
        ax1.axvline(x=crossed_lower_bound_at, lw=setlw*5, color='black') # plot the x value where the m value first enters the 95 percentile
        ax1.text(x=crossed_lower_bound_at, y=0.4, ha='left', va='center', s=crossed_text, bbox=bbox, fontsize=text_size)

        # plot normal distribution next to m plot
        ax2.plot(normpdf, pdf_space, lw=setlw*5, color='red', label=norm_label) # plot the normal distribution of the m values along the y axis
        
        ax1.legend(facecolor=fc, framealpha=box_alpha, edgecolor=ec)#, fontsize=text_size)
        ax2.legend(facecolor=fc, framealpha=box_alpha, edgecolor=ec)#, fontsize=text_size)
        ax0.set_xlabel("inversion cycle", fontsize=text_size)
        ax0.set_ylabel("unique interactions", fontsize=text_size)
        ax1.set_xlabel("inversion cycle", fontsize=text_size)
        ax1.set_ylabel(r"$m$", fontsize=text_size)
        fig.suptitle(plot_title, fontsize=plot_title_size)
        ax0.set_title(subplot0_title, fontsize=subplot_title_size)
        ax1.set_title(subplot1_title, fontsize=subplot_title_size)
        #ax3.set_title(subplot3_title, fontsize=subplot_title_size)

        # create the diagram directory if it does not exist yet

        output_dir=self.get_LiSC_path('diagrams/')
        output_name='inversion_sim_a{a}_b{b}'.format(a=self.genesA, b=self.genesB)
            
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        if not os.path.exists(output_dir+'/pdf'):
            os.makedirs(output_dir+'/pdf')
        if not os.path.exists(output_dir+'/png'):
            os.makedirs(output_dir+'/png')
        
        # save this as a pdf and png
        plt.savefig(output_dir+'pdf/'+output_name+'.pdf')
        plt.savefig(output_dir+'png/'+output_name+'.png')

        if not yaml:
            return
        
        # save the trace as a yaml file
        import yaml
        if not os.path.exists(output_dir+'/yaml'):
            os.makedirs(output_dir+'/yaml')
        with open(output_dir+'yaml/'+output_name+'.yaml', 'w') as f:
            yaml.dump(self.trace, f)

    def log(self, elapsed='-1'):
        output_dir=self.get_LiSC_path('log/')
        log_file='log.csv'

        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        newfile=not os.path.exists(output_dir+log_file)
        mode='w' if newfile else 'a'
        
        header='timestamp;|A|;|B|;cycles;t50;AB_convergence;first_95_m;m_sigma;m_mu;Delta_t\n'
        format_string='{ts};{A};{B};{c};{t50};{ABconv};{m95:.3f};{sig:.3f};{mu:.3f};{dt}\n'

        from datetime import datetime as dt
        
        with open(output_dir+log_file, mode) as f:
            if newfile:
                f.write(header)
            f.write(format_string.format(ts=str(dt.now()), A=self.genesA, B=self.genesB, c= self.cycle, t50=self.t50, ABconv=self.AB_convergence, m95=self.first_95_m, sig=self.m_sigma, mu=self.m_mu, dt=elapsed))

    """
    get the proper path if the program is running on LiSC
    """
    def get_LiSC_path(self, output_dir):
        return (remote_path_base if os.path.exists(remote_path_base) else '')+output_dir
        if os.path.exists(remote_path_base):
            return remote_path_base+output_dir
        return output_dir
        
def print_usage():
    message="""
    usage: inversion_sim.py <Asize> <Bsize>

    - Asize: number of genes in group A (int)
    - Bsize: number of genes in group B (int) 
    """
    print(message)
            
def main():
    
    iterations = 100#000
    
    size_pairs=[(5, 5), (5, 10), (10, 100), (500, 400), (1000, 1000), (500, 1000)] # pairs of A and B sizes to simulate

    # handle command line arguments
    if len(sys.argv) != 3:
        print_usage()
        return

    Asize=-1
    Bsize=-1
    try:
        Asize=int(sys.argv[1])
        Bsize=int(sys.argv[2])
    except ValueError:
        print_usage()

    # start a timer
    start=time.time()

    print("\ncreating chromosome...")
    chrom=Chrom(Asize+Bsize, Asize, Bsize)
    #chroms=[Chrom(10000000, pair[0], pair[1]) for pair in size_pairs] # create chromosomes
    print("\nrunning simulation...")
    chrom.simulation_cycle(until_converged=True)
    #for chrom in chroms: # run all simulations
    #    chrom.simulation_cycle(iterations=iterations//max([pair[0]+pair[1] for pair in size_pairs]*(chrom.genesA+chrom.genesB))*(chrom.genesA+chrom.genesB))
    print("\nplotting results...")
    chrom.plot_results()
    #for chrom in chroms: # plot all results
    #    print("plotting |A|={A:4d}, |B|={B:4d}".format(A=chrom.genesA, B=chrom.genesB))
    #    chrom.plot_results()

    end=time.time()
    elapsed=end-start
    elapsed_string="{minutes:02d}:{seconds:02d}".format(minutes=int(elapsed//60), seconds=int(elapsed%60))
    print("\nelapsed time: "+elapsed_string)
    chrom.log(elapsed_string)
    
if __name__ == "__main__":
    main()
