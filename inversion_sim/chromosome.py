"""
This program calculates the time it takes for any gene to interact with another
after a fusion-with mixing event.
"""

import random
import numpy as np

class Chrom():
    def __init__(self, length, gene_quantityA, gene_quantityB, level_of_convergence=1, window_size=1, until_converged=False):
        # length is the chromosome length
        # The intention is to eventually model varying regions of gene density,
        #    so don't delete this yet.
        self.length = length
        self.genesA = gene_quantityA
        self.genesB = gene_quantityB
        self.size=self.genesA+self.genesB

        self.gene_list = ["A."+str(i+1) for i in range(self.genesA)] + ["B."+str(i+1) for i in range(self.genesB)]
        # seen is a list of genes that have already interacted with each other, value is cycle
        self.seen = {}
        self.window_size=window_size
        self.t50=-1
        self.AB_convergence=-1
        self.converged_AtoB=0
        self.converged_BtoA=0
        # level_of_convergence specififies the fraction of possible gene interactions have
        #   to have been seen before ending the simulation if until_converged is set to True.
        # The value ranges from 0 to 1.
        self.level_of_convergence=level_of_convergence
        self.until_converged=until_converged

        # set the rate of data sampling (each cycle for small chromosomes, all the way to 1% of cycles for large ones)
        self.sample_frequency = 1 if self.size <= 50 else 0.5 if self.size <= 100 else 0.1 if self.size <= 500 else 0.01
        self.sample_rate = int( 1 / self.sample_frequency)

        # this cumulatively tracks how many new interactions were added at each cycle
        self.trace      = {k:[] for k in self.gene_list}
        self.trace_AtoB = {k:[] for k in self.gene_list if k.startswith("A")}
        self.trace_BtoA = {k:[] for k in self.gene_list if k.startswith("B")}
        self.trace_m={}
        self.cycle = 0
        # [PERFORMANCE] Parallelizing these two functions could save some time.
        self.update_seen()
        self.calculate_m()
        self.first_95_m=-1
        self.m_sigma=-1
        self.m_mu=-1
    
    def __str__(self):
        format_string="cycle: {cycle:10d}; last inversion: {start} {istart:5d} {inverted} {iend:<5d} {end}; m: {m:1.3f}"
        AB_string=self.get_AB_string()
        i0=self.last_i0
        i1=self.last_i1
        m=self.trace_m[self.cycle]
        ret=format_string.format(cycle=self.cycle, start=AB_string[0:i0], istart=i0, inverted=AB_string[i0:i1], iend=i1-1, end=AB_string[i1:len(AB_string)], m=m)
        return ret
    
    def simulation_cycle(self, iterations = 0):
        """
        This function runs the simulation for a given number of iterations, or until done.
        """
        # raise an error if we don't know how to run this method
        if iterations == 0 and self.until_converged == False:
            raise ValueError("Please specify either iterations or until_converged.")
        if self.until_converged:
            converging_at=self.calculate_convergence()
            AB_have_converged=False
            while len(self.seen)/converging_at < self.level_of_convergence: 
                self.shuffle()

                if self.t50 < 0 and len(self.seen) >= converging_at/2:
                    self.t50=self.cycle
                    print("reached t50 at "+str(self.t50))
                if self.t50 >= 0 and not AB_have_converged:
                    AB_have_converged=self.converged_AtoB >= self.genesA and self.converged_BtoA >= self.genesB
                    if AB_have_converged:
                        self.AB_convergence=self.cycle
                        print("A-B/B-A converged at "+str(self.AB_convergence))
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

        # calculate all the values        
        cycles=[x for x in self.trace_m.keys()]
        m_values=[y for y in self.trace_m.values()]
        burn_in=0.25
        start_norm_at=int((self.cycle+1)*burn_in)//1
        burnt_m_values=np.array(m_values[start_norm_at:])
        mu=np.mean(burnt_m_values)
        var=np.var(burnt_m_values)
        sigma=var**0.5
        self.m_sigma=sigma
        self.m_mu=mu
        upper_bound=mu+1.96*sigma
        lower_bound=mu-1.96*sigma
        crossed_lower_bound_at=0
        for k in self.trace_m.items():
            if k[1] >= lower_bound:
                crossed_lower_bound_at=k[0]
                break
        self.first_95_m=crossed_lower_bound_at   

    def calculate_convergence(self):
        """
        calculates the number of possible unique gene interactions
        """
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
        # [PERFORMANCE] Maybe parallelizing update_seen() and calculate_m() here would save some time.
        self.update_seen()
        self.cycle += 1
        self.calculate_m()

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
            #print(range(i+1, min(len(self.gene_list)-1, i+self.window_size+1)))
            for l in range(i+1, min(len(self.gene_list), i+self.window_size+1)):
                #print("in loop")
                this_edge = tuple(sorted([self.gene_list[i], self.gene_list[l]]))
                if this_edge not in self.seen:
                    self.seen[this_edge] = self.cycle
                    # update the trace structure
                    for j in [0,1]:
                        other = 1 if j == 0 else 0
                        self.trace[this_edge[j]][-1] += 1
                        if this_edge[j].startswith("A") and this_edge[other].startswith("B"):
                            self.trace_AtoB[this_edge[j]][-1] += 1
                            # check whether the latest count of cross-group gene interactions for this gene is equal
                            #   to the number of genes in the opposite group, increase counter for converged A-to-B genes
                            # subtract 1 from the needed interactions if this gene is A.1 (telomeres stay intact and
                            #   cannot interact with the other telomere)
                            if self.trace_AtoB[this_edge[j]][-1] == self.genesB-(1 if this_edge[j].split('.')[1] == '1' else 0):
                                self.converged_AtoB+=1
                        if this_edge[j].startswith("B") and this_edge[other].startswith("A"):
                            self.trace_BtoA[this_edge[j]][-1] += 1
                            # same as above, but for B-to-A
                            if self.trace_BtoA[this_edge[j]][-1] == self.genesA-(1 if this_edge[j].split('.')[1] == str(self.genesB) else 0):
                                self.converged_BtoA+=1
                
    def get_AB_string(self):
        return ''.join([gene[0] for gene in self.gene_list])

    def calculate_m(self):
        """
        calculate m of the current gene sequence
        """
        # [PERFORMANCE] This looks like it runs a lot of loops in the background. Maybe we could save some runtime by running one loop
        #   and doing the actions manually?
        sequence=self.get_AB_string()
        substrings = [sequence[i:i+2] for i in range(len(sequence)-1)]
        A = sequence.count('A') # These two lines seem reduntant, since we have the A and B counts stored in the chromosome.
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
