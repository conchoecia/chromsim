"""
This file contains the class for a chromosome with the ability to simulate inversions.
"""

import random
import numpy as np
import datetime as dt

def check_AB_pair(AB0, AB1):
    return 0 if AB0[0] == AB1[0] else 1

class Chrom():
    def __init__(self, Asize, Bsize, length=0, level_of_convergence=1, window_size=1, inversion_cuts=[], timestamp=None):

        # set parameters
        self.Asize=Asize
        self.Bsize=Bsize
        self.size=self.Asize+self.Bsize
        self.length=length
        self.level_of_convergence=level_of_convergence
        self.window_size=window_size

        # set the initial list of cuts
        self.inversion_cuts=inversion_cuts if inversion_cuts else []
        
        # set constants based on parameters
        self.m_const=2*self.Asize*self.Bsize/(self.Asize+self.Bsize-1)
        self.converging_at=self._calculate_convergence()

        # set the cycle number to 0
        self.cycle=0

        # initialize gene list
        self.gene_list = ["A."+str(i+1) for i in range(self.Asize)] + ["B."+str(i+1) for i in range(self.Bsize)]
        self.original_gene_list=self.gene_list.copy()
        self.t50_gene_list=None

        # initialize the points of (half-)convergence and entropy
        self.t50=-1
        self.t100=-1
        self.AB_convergence=-1 # t100 but for inter-group convergence only (don't care about intra-group convergence)
        self.tS=-1
        self.m_sigma=-1
        self.m_mu=-1
        
        # initialize the number of converged genes in both groups (i.e. the number of genes in group A that has seen all group B genes, and vice versa
        self.converged_AtoB=0 # once these two numbers reach Asize and Bsize, respectively,
        self.converged_BtoA=0 # AB_convergence is achieved
        
        # initialize the traces that cumulatively track how many new interactions are added after each cycle
        self.trace={k:{0:(1 if k == 'A.1' or k == 'B.'+str(self.Bsize+1) else 2)} for k in self.gene_list}
        self.last_trace_values={k:self.trace[k][0] for k in self.trace}
        self.trace_mean={0:np.mean(list(self.last_trace_values.values()))}

        self.trace_AtoB={k:{0:(1 if k == 'A.'+str(self.Asize) else 0)} for k in self.gene_list if k.startswith("A")}
        self.last_trace_AtoB_values={k:self.trace_AtoB[k][0] for k in self.trace_AtoB}
        self.trace_AtoB_mean={0:np.mean(list(self.last_trace_AtoB_values.values()))}

        self.trace_BtoA={k:{0:(1 if k == 'B.1' else 0)} for k in self.gene_list if k.startswith("B")}
        self.last_trace_BtoA_values={k:self.trace_BtoA[k][0] for k in self.trace_BtoA}
        self.trace_BtoA_mean={0:np.mean(list(self.last_trace_BtoA_values.values()))}

        # initialize the entropy value trace
        self.trace_m={0:0}

        # initialize dictionary of new seen interactions
        self.seen = {} # seen is a list of genes that have already interacted with each other, value is the cycle of the first interaction between the two genes that are the key
        self._update_seen(0, len(self.gene_list)-1, trace=False)

        # give this chromosome a timestamp
        self.timestamp=timestamp if timestamp else dt.datetime.now()    

    def _print_progress(self, n):
        """
        print a line showing the progess of the simulation

        s is appended to the end of the string
        """
        
        print("\rcycle {cycle:15d} {progress:.2f}%".format(cycle=self.cycle, progress=(self.cycle/n if n >= 0 else len(self.seen)/self.converging_at)*100), end="")

    def _quit_condition(self, n, m):
        """
        return True or False depending on whether the condition to end the simulation is met (convergence/cycles/m)
        """

        return self.trace_m[self.cycle] >= m if m >= 0 else self.cycle >= n if n >= 0 else len(self.seen) >= self.level_of_convergence*self.converging_at

    def run(self, n=-1, m=-1, show_output=False, trace=True):
        """
        run the simulation for n iterations, until a given m is reached, or until convergence if n,m < 0

        print progress report to the console if show_output is True
        """

        if n >= 0 and m >= 0:
            raise Exception('n and m cannot both be specified for a simulation run.')
        
        # run the simulation for the specified number of iterations/until convergence
        while not self._quit_condition(n, m):
            if show_output and self.cycle%100 == 0:
                self._print_progress(n)
                
            self._shuffle(trace)
            
            if self.t50 < 0 and len(self.seen) >= self.converging_at/2:
                self.t50=self.cycle
                self.t50_gene_list=self.gene_list.copy()
            if self.t100 < 0 and len(self.seen) >= self.converging_at:
                self.t100=self.cycle
                self.t100_gene_list=self.gene_list.copy()
            if self.t50 >= 0 and self.AB_convergence < 0:
                if self.converged_AtoB >= self.Asize and self.converged_BtoA >= self.Bsize:
                    self.AB_convergence=self.cycle

        if show_output:
            # show the progress bar with 100% completion
            self._print_progress(n)
        
        self._calculate()
        
    def _calculate(self):
        """
        calculate all the values based on the simulation results       
        """

        # get the cycles and m_values as lists
        cycles=[x for x in self.trace_m.keys()]
        m_values=[y for y in self.trace_m.values()]

        # set the burn-in to 25% of the data
        burn_in=0.25
        start_norm_at=int((self.cycle+1)*burn_in)//1
        burnt_m_values=np.array(m_values[start_norm_at:])

        # calculate μ and σ
        mu=np.mean(burnt_m_values)
        var=np.var(burnt_m_values)
        sigma=var**0.5 # σ is the square root of the variance
        self.m_sigma=sigma
        self.m_mu=mu

        # set the bounds of the 5th and 95th percentiles (values in this range are considered entropic)
        upper_bound=mu+1.96*sigma
        lower_bound=mu-1.96*sigma

        # find the first entropic value
        tS=0
        for k in self.trace_m.items():
            if k[1] >= lower_bound:
                tS=k[0]
                break
        self.tS=tS   
        
    def _calculate_convergence(self):
        """
        calculates the number of possible unique gene interactions
        """

        conv=0
        for i in range(self.Asize+self.Bsize):
            conv+=i
        return conv-1

    def _pick_cut(self):
        """
        return a tuple of indices that represent the breakpoints of the next inversion
        """

        if self.cycle < len(self.inversion_cuts):
            return self.inversion_cuts[self.cycle]

        valid=False
        while not valid:
            cut=(random.randint(1, len(self.gene_list)-1), random.randint(1, len(self.gene_list)-1))
            valid=self._validate_cut(cut)
        return cut

    def _validate_cut(self, cut):
        """
        check whether the two breakpoints chosen should be rejected or not

        no rules are imposed currently, returns True
        """
        
        return True

        # discard the cut if a random number is above 1/(i0-i1) => longer inversions are less likely
        import random as r
        distance=abs(cut[1]-cut[0])
        if distance == 0:
            return False
        cutoff=1/distance
        rand=r.random()
        return rand <= cutoff
            
    def _shuffle(self, trace):
        """
        invert a random chromosome section
        """

        self._invert()
        self._update_cycle(self.inversion_cuts[-1], trace)

    def _invert(self):
        """
        pick cut and invert the chromosome
        """
        
        cut=sorted(self._pick_cut())
        self.inversion_cuts.append(cut)

        i0 = cut[0]
        i1 = cut[1]
        self.gene_list[i0:i1] = self.gene_list[i0:i1][::-1]
        
    def _get_window(self, i):
        """
        get the window around an index, based upon window_size
        """

        # the max and min operations make sure that there is no out-of-bounds action going on
        return max(i-self.window_size, 0), min(i+self.window_size, len(self.gene_list)-1)
        
    def _update_cycle(self, cut, trace):
        """
        prepare and execute updates of seen and m (only in the window around the breakpoints of cut)
        """

        i0 = cut[0]
        i1 = cut[1]

        # increment cycle
        self.cycle+=1

        # update seen around i0 and i1
        start, end=self._get_window(i0)
        self._update_seen(start, end, trace)
        start, end=self._get_window(i1)
        self._update_seen(start, end, trace)

        # update m around i0 and i1
        self._update_m(i0, i1)
        
    def _update_seen(self, start, end, trace):
        """
        update the seen list and trace structures

        start/end: indices demarking the window where the previous changes are relevant
        trace: boolean value specifying whether or not to save data to the trace structures
        """
        
        # iterate over all genes between start and end
        for gene_index in range(start, end):
            gene0=self.gene_list[gene_index]
            
            # iterate over all following genes to compare to the current one
            for comparing_gene_index in range(gene_index+1, min(end+1, gene_index+self.window_size+1)):
                # group the two genes in a tuple
                gene1=self.gene_list[comparing_gene_index]
                pair=tuple(sorted([gene0, gene1]))
                
                # update data structures if the pairing is new
                if pair not in self.seen:
                    # udpate seen with the current cycle
                    self.seen[pair]=self.cycle

                    # skip the following steps if trace is False
                    if not trace:
                        continue
                    # update trace structure for both genes in the tuple
                    for this_index in [0,1]:
                        other_index=1 if this_index == 0 else 0

                        this=pair[this_index]
                        other=pair[other_index]
                        
                        # update the current gene's trace with the current cycle and increment the last entry by 1 => (cycle, number_of_interactions)
                        if self.cycle > 0:
                            previous_value=self.last_trace_values[this]
                            self.trace[this][self.cycle]=previous_value+1
                            self.last_trace_values[this]=previous_value+1
                            self.trace_mean[self.cycle]=np.mean([v for v in self.last_trace_values.values()])

                        # update inter-group traces with a new interaction
                        if this.startswith("A") and other.startswith("B"):
                            if self.cycle > 0:
                                previous_value=self.last_trace_AtoB_values[this]
                                self.trace_AtoB[this][self.cycle]=previous_value+1
                                self.last_trace_AtoB_values[this]=previous_value+1
                                self.trace_AtoB_mean[self.cycle]=np.mean([v for v in self.last_trace_AtoB_values.values()])
                            # check whether the latest count of cross-group gene interactions for this gene is equal
                            #   to the number of genes in the opposite group, increase counter for converged A-to-B genes
                            # subtract 1 from the needed interactions if this gene is A.1 (telomeres stay intact and
                            #   cannot interact with the other telomere)
                            if self.trace_AtoB[pair[this_index]][self.cycle] == self.Bsize-(1 if this.split('.')[1] == '1' else 0):
                                self.converged_AtoB+=1
                        if this.startswith("B") and other.startswith("A"):
                            if self.cycle > 0:
                                previous_value=self.last_trace_BtoA_values[this]
                                self.trace_BtoA[this][self.cycle]=previous_value+1
                                self.last_trace_BtoA_values[this]=previous_value+1
                                self.trace_BtoA_mean[self.cycle]=np.mean([v for v in self.last_trace_BtoA_values.values()])
                            # same as above, but for B-to-A
                            if self.trace_BtoA[this][self.cycle] == self.Asize-(1 if this.split('.')[1] == str(self.Bsize) else 0):
                                self.converged_BtoA+=1
                
    def _update_m(self, i0, i1):
        """
        update m based on changes around the breakpoints
        """

        # get the previous m value
        old_m=self.trace_m[self.cycle-1]

        # do nothing if the inversion was of length 0
        if i0==i1:
            self.trace_m[self.cycle]=old_m
            return

        # calculate the number of transition before the last inversion
        old_transitions=old_m*self.m_const+1

        # calculate the difference in the number of transtions between the old and new pairings
        old_pair0=(self.gene_list[i0-1], self.gene_list[i1-1])
        old_pair1=(self.gene_list[i0], self.gene_list[i1])
        new_pair0=(self.gene_list[i0-1:i0+1])
        new_pair1=(self.gene_list[i1-1:i1+1])
        delta_transitions=check_AB_pair(new_pair0[0], new_pair0[1])+check_AB_pair(new_pair1[0], new_pair1[1])-check_AB_pair(old_pair0[0], old_pair0[1])-check_AB_pair(old_pair1[0], old_pair1[1])

        # caluculate the new number of transitions
        new_transitions=old_transitions+delta_transitions

        # calculate the new m value and append to the data structure
        new_m=(new_transitions-1)/self.m_const
        self.trace_m[self.cycle]=new_m
