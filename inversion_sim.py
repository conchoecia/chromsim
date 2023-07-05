#!/usr/bin/env python

"""
This program calculates the time it takes for any gene to interact with another
after a fusion-with mixing event.
"""

import random
import sys


class Chrom():
    def __init__(self, size, gene_quantityA, gene_quantityB):
        self.length = size
        self.genesA = gene_quantityA
        self.genesB = gene_quantityB
        self.gene_list = ["A."+str(i+1) for i in range(self.genesA)] + ["B."+str(i+1) for i in range(self.genesB)]
        # seen is a list of genes that have already interacted with each other, value is cycle
        self.seen = {}
        self.sample_frequency = .01 # this is the percent of the data that we sample, range [0-1]. 0.01 is a good rate for samples above 100k iterations
        self.sample_rate = int( 1 / self.sample_frequency)
        print("sample rate: every ", self.sample_rate, " iterations sampled")
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
                print("\r{0:15d}  {1:.2f}%  ".format(i, (i/iterations)*100), end="")
            self.shuffle()
        print("{0:15d}  100.00%  ".format(i+1))
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
        
        # initialize the matplotlib plot
        import matplotlib.pyplot as plt

        #print([k for k in  self.trace_AtoB])
        #print("A.1: ", self.trace_AtoB["A.1"])
        #print([k for k in  self.trace_BtoA], self.trace_BtoA)

        setlw = 0.2
        seta  = 0.1
        A_alpha   = max(1/self.genesA, 0.05)
        B_alpha   = max(1/self.genesB, 0.05)
        all_alpha = max(1/(self.genesA + self.genesB), 0.05)

        # set up the panel
        for k in self.trace:
            plt.plot([x*self.sample_rate for x in range(len(self.trace[k]))],
                     self.trace[k], color='black', lw = setlw, alpha=all_alpha)
        plt.plot([x*self.sample_rate for x in range(len(self.trace[k]))],
                 self._median_of_trace(self.trace), color='black', lw = setlw*10, alpha=0.75)
 
        # now plot the A-to-B and B-to-A traces
        for k in self.trace_AtoB:
            plt.plot([x*self.sample_rate for x in range(len(self.trace_AtoB[k]))],
                     self.trace_AtoB[k], color='blue', lw = setlw, alpha=A_alpha)
        plt.plot([x*self.sample_rate for x in range(len(self.trace_AtoB[k]))],
                 self._median_of_trace(self.trace_AtoB), color='blue', lw = setlw*10, alpha=0.75)


        for k in self.trace_BtoA:
            plt.plot([x*self.sample_rate for x in range(len(self.trace_BtoA[k]))],
                     self.trace_BtoA[k], color='red', lw = setlw, alpha=B_alpha)
        plt.plot([x*self.sample_rate for x in range(len(self.trace_BtoA[k]))],
                 self._median_of_trace(self.trace_BtoA), color='red', lw = setlw*10, alpha=0.75)

 
        plt.xlabel("inversion cycle")
        plt.ylabel("Unique interactions")
        # save this as a pdf
        plt.savefig("inversion_sim.pdf")

        # save the trace as a yaml file
        import yaml
        with open("inversion_sim.yaml", "w") as f:
            yaml.dump(self.trace, f)

def main():
    iterations = 100000
    chrom = Chrom(10000000, 500, 400)
    chrom.simulation_cycle(iterations = iterations)
    chrom.plot_results()

if __name__ == "__main__":
    main()
