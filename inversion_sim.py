#!/usr/bin/env python

"""
This program calculates the time it takes for any gene to interact with another
after a fusion-with mixing event.
"""

import random
import sys


class Chrom():
    def __init__(self, size, gene_quantityA, gene_quantityB):
        self.length = size # is this even relevant when gene_quantityA and gene_quantityB are specified?
        self.genesA = gene_quantityA
        self.genesB = gene_quantityB
        self.gene_list = ["A."+str(i+1) for i in range(self.genesA)] + ["B."+str(i+1) for i in range(self.genesB)]
        # seen is a list of genes that have already interacted with each other, value is cycle
        self.seen = {}
        self.sample_frequency = .01 # this is the percent of the data that we sample, range [0-1]. 0.01 is a good rate for samples above 100k iterations
        self.sample_rate = int( 1 / self.sample_frequency)
        #print("sample rate: every ", self.sample_rate, " iterations sampled")
        # this tracks how many new interactions were added at each cycle
        self.trace      = {k:[] for k in self.gene_list}
        self.trace_AtoB = {k:[] for k in self.gene_list if k.startswith("A")}
        self.trace_BtoA = {k:[] for k in self.gene_list if k.startswith("B")}
        self.trace_m={}
        self.cycle = 0
        self.update_seen()
        self.calculate_m()
    
    def simulation_cycle(self, iterations = 0, until_converged = False):
        """
        This function runs the simulation for a given number of iterations, or until done.
        """
        # raise an error if we don't know how to run this method
        if iterations == 0 and until_converged == False:
            raise ValueError("Please specify either iterations or until_converged.")
        # until_converged hasn't been implemented yet, so tell the user that this isn't possible now
        if until_converged:
            raise NotImplementedError("This function is not yet implemented.")
        # run the simulation for the specified number of iterations
        for i in range(iterations):
            # make a progress bar that stays on the same line that prints every 100 cycles
            if i % 100 == 0:
                # use a carriage return to stay on the same line
                # use format string to make the number occupy 10 spaces.
                # separate with a space character, then plot the percentage.
                print("\r|A|={Asize:5d}, |B|={Bsize:5d}: {cycle:15d} {progress:.2f}%  ".format(Asize=self.genesA, Bsize=self.genesB, cycle=i, progress=(i/iterations)*100), end="")
            self.shuffle()
        print("{cycle:15d}  100.00%  ".format(cycle=i+1))
        print("")
 
    def shuffle(self):
        # Randomly pick two indices in the list.
        # Start at one and end at len-1 to not destroy telomeres
        i1 = random.randint(1, len(self.gene_list)-1)
        i2 = random.randint(1, len(self.gene_list)-1)
        sortedi = sorted([i1, i2]) 
        i1 = sortedi[0]
        i2 = sortedi[1]
        self.gene_list[i1:i2] = self.gene_list[i1:i2][::-1]
        self.update_seen()
        self.calculate_m()
        self.cycle += 1
    
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

    """
    calculate m of the current gene sequence
    """
    def calculate_m(self):
        sequence=''.join([gene[0] for gene in self.gene_list])
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

    def plot_results(self):
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
        start_norm_at=math.floor(self.cycle*burn_in)
        burnt_m_values=np.array(m_values[start_norm_at:])
        mu=np.mean(burnt_m_values)
        var=np.var(burnt_m_values)
        sigma=math.sqrt(var)
        pdf_space=np.linspace(0, 2, 100)
        normpdf=stats.norm.pdf(pdf_space, mu, sigma)
        scaled_normpdf=normpdf/max(normpdf)*self.cycle/10
        upper_bound=mu+1.96*sigma
        lower_bound=mu-1.96*sigma
        crossed_lower_bound_at=0
        for k in self.trace_m.items():
            if k[1] >= lower_bound:
                crossed_lower_bound_at=k[0]
                break
        cycle_limit=start_norm_at*2
        sample_limit=cycle_limit//self.sample_rate

        # initialize the matplotlib plot
        import matplotlib.pyplot as plt
        plt.style.use('bmh')

        #print([k for k in  self.trace_AtoB])
        #print("A.1: ", self.trace_AtoB["A.1"])
        #print([k for k in  self.trace_BtoA], self.trace_BtoA)

        setlw = 0.2
        seta  = 0.1
        A_alpha   = max(1/self.genesA, 0.05)
        B_alpha   = max(1/self.genesB, 0.05)
        all_alpha = max(1/(self.genesA + self.genesB), 0.05)
        figsize=(10, 15)
        bbox=dict(facecolor = 'white', alpha=0.5, boxstyle='round, pad=0.5', edgecolor='gray')
        plot_title_size=20
        subplot_title_size=15

        plot_title=r"$|A|={Asize}, |B|={Bsize}$, {cycles} inversion cycles".format(Asize=self.genesA, Bsize=self.genesB, cycles=self.cycle)
        subplot0_title=r"number of unique interactions after $n$ inversion cycles"
        subplot1_title=r"$m$ after $n$ inversion cycles"
        
        fig, axes=plt.subplots(2, 1, figsize=figsize)
        
        # set up the panel
        axes[0].set_xlim(0, self.cycle)
        for k in self.trace:
            axes[0].plot([x*self.sample_rate for x in range(len(self.trace[k]))],
                     self.trace[k], color='black', lw = setlw, alpha=all_alpha)
        axes[0].plot([x*self.sample_rate for x in range(len(self.trace[k]))],
                 self._median_of_trace(self.trace), color='black', lw = setlw*10, alpha=0.75)
 
        # now plot the A-to-B and B-to-A traces
        for k in self.trace_AtoB:
            axes[0].plot([x*self.sample_rate for x in range(len(self.trace_AtoB[k]))],
                     self.trace_AtoB[k], color='blue', lw = setlw, alpha=A_alpha)
        axes[0].plot([x*self.sample_rate for x in range(len(self.trace_AtoB[k]))],
                 self._median_of_trace(self.trace_AtoB), color='blue', lw = setlw*10, alpha=0.75)


        for k in self.trace_BtoA:
            axes[0].plot([x*self.sample_rate for x in range(len(self.trace_BtoA[k]))],
                     self.trace_BtoA[k], color='red', lw = setlw, alpha=B_alpha)
        axes[0].plot([x*self.sample_rate for x in range(len(self.trace_BtoA[k]))],
                 self._median_of_trace(self.trace_BtoA), color='red', lw = setlw*10, alpha=0.75)
        
        #plot m values on the second subplot
        axes[1].set_xlim([0, cycle_limit])
        axes[1].plot(cycles[:cycle_limit], m_values[:cycle_limit], lw=setlw*2, color='blue', label=r"$m$")

        # plot 95 percentile of m value normal distribution
        crossed_text="first 95 percentile value:\n{cross} cycles\n({perc:.2f}% of cycles)".format(cross=crossed_lower_bound_at, perc=crossed_lower_bound_at/self.cycle*100)
        burn_in_text="burn-in:\n{cycle} cycles\n({perc:.0f}% of cycles)".format(cycle=start_norm_at, perc=burn_in*100)
        norm_label=r"normal distribution of $m$" "\n" "(excluding the first {perc}% of cycles)".format(perc=burn_in*100)
        axes[1].axhline(y=upper_bound, color='red', lw=setlw*5, ls=':') # plot upper bound of the 95 percentile
        axes[1].axhline(y=lower_bound, color='red', lw=setlw*5, ls=':') # plot lower bound of the 95 percentile
        axes[1].text(y=upper_bound, x=cycle_limit, ha='right', va='bottom', s=r"$\mu+1.96\cdot\sigma$", bbox=bbox)
        axes[1].text(y=lower_bound, x=cycle_limit, ha='right', va='top', s=r"$\mu-1.96\cdot\sigma$", bbox=bbox)
        axes[1].fill_between(cycles, lower_bound, upper_bound, color='red', alpha=0.1) # shade area between the bounds of the 95 percentile
        #axes[1].axvline(x=start_norm_at, lw=setlw*5, color='black') # plot the x value of the burn-in
        axes[1].axvline(x=crossed_lower_bound_at, lw=setlw*5, color='black') # plot the x value where the m value first enters the 95 percentile
        #axes[1].text(x=start_norm_at, y=0.2, ha='left', va='center', s=burn_in_text, bbox=bbox)
        axes[1].text(x=crossed_lower_bound_at, y=0.4, ha='left', va='center', s=crossed_text, bbox=bbox)
        axes[1].plot(scaled_normpdf, pdf_space, lw=setlw*5, color='red', label=norm_label) # plot the normal distribution of the m values along the y axis

        axes[1].legend(facecolor='white', framealpha=0.5, edgecolor='gray')
        plt.xlabel("inversion cycle")
        axes[0].set_ylabel("unique interactions")
        axes[1].set_ylabel(r"$m$")
        fig.suptitle(plot_title, fontsize=plot_title_size)
        axes[0].set_title(subplot0_title)
        axes[1].set_title(subplot1_title)

        output_dir='diagrams/'
        output_name='inversion_sim_a{a}_b{b}'.format(a=self.genesA, b=self.genesB)

        # create the diagram directory if it does not exist yet
        import os
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        if not os.path.exists(output_dir+'/pdf'):
            os.makedirs(output_dir+'/pdf')
        if not os.path.exists(output_dir+'/png'):
            os.makedirs(output_dir+'/png')
        if not os.path.exists(output_dir+'/yaml'):
            os.makedirs(output_dir+'/yaml')
        
        # save this as a pdf and png
        plt.savefig(output_dir+'pdf/'+output_name+'.pdf')
        plt.savefig(output_dir+'png/'+output_name+'.png')

        # save the trace as a yaml file
        import yaml
        with open(output_dir+'yaml/'+output_name+'.yaml', "w") as f:
            yaml.dump(self.trace, f)

def main():
    iterations = 100000
    #chrom = Chrom(10000000, 500, 400)
    
    size_pairs=[(5, 5), (5, 10), (10, 100), (500, 400), (1000, 1000), (500, 1000)] # pairs of A and B sizes to simulate
    print("creating chromosomes...")
    chroms=[Chrom(10000000, pair[0], pair[1]) for pair in size_pairs] # create chromosomes
    print("running simulations...")
    for chrom in chroms: # run all simulations
        chrom.simulation_cycle(iterations=iterations//max([pair[0]+pair[1] for pair in size_pairs]*(chrom.genesA+chrom.genesB))*(chrom.genesA+chrom.genesB))
    print("plotting results...")
    for chrom in chroms: # plot all results
        chrom.plot_results()


if __name__ == "__main__":
    main()
