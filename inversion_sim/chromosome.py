"""
This file contains the class for a chromosome with the ability to simulate inversions.
"""

import random
import numpy as np

class Chrom():
    def __init__(self, length, Asize, Bsize, level_of_convergence=1, window_size=1, translocations_per_cycle=0, cuts=[]):

        # set parameters
        self.length=length
        self.genesA=Asize
        self.genesB=Bsize
        self.size=self.genesA+self.genesB
        self.level_of_convergence=level_of_convergence
        self.translocations_per_cycle=translocations_per_cycle
        self.window_size=window_size
        self.inversion_cuts=cuts if cuts else []
        
        # set constants based on parameters
        self.m_const=2*self.genesA*self.genesB/(self.genesA+self.genesB-1)
        self.converging_at=self._calculate_convergence()

        # set the cycle number to 0
        self.cycle=0

        # initialize gene list
        self.gene_list = ["A."+str(i+1) for i in range(self.genesA)] + ["B."+str(i+1) for i in range(self.genesB)]
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
        self.trace      = {k:[(0, 0)] for k in self.gene_list}
        self.trace_AtoB = {k:[(0, 0)] for k in self.gene_list if k.startswith("A")}
        self.trace_BtoA = {k:[(0, 0)] for k in self.gene_list if k.startswith("B")}

        # initialize the entropy value trace
        self.trace_m={0:0}

        # initialize dictionary of new seen interactions
        self.seen = {} # seen is a list of genes that have already interacted with each other, value is the cycle of the first interaction between the two genes that are the key
        self._update_seen(0, len(self.gene_list)-1)

    def _print_progress(self, progress, s=""):
        """
        print a line showing the progess of the simulation

        s is appended to the end of the string
        """
        print("\rcycle {cycle:15d} {progress:.2f}%".format(cycle=self.cycle, progress=progress), end="")
        
    def run(self, n=-1, show_output=True):
        """
        run the simulation for n iterations, or until convergence if n < 0

        print progress report to the console if show_output is True
        """

        # run n iterations if n >= 0
        if n >= 0:
            # run the simulation for the specified number of iterations
            while self.cycle < n:
                if show_output and self.cycle%100 == 0:
                    self._print_progress(self.cycle/n*100)
                self._shuffle()
            if show_output:
                # show the progress bar with 100% completion
                self._print_progress(100)
        # run until convergence otherwise
        else:
            AB_have_converged=False
            while len(self.seen)/self.converging_at < self.level_of_convergence:
                if show_output and self.cycle%100 == 0:
                    self._print_progress(len(self.seen)/self.converging_at*100)
                self._shuffle()
                
                if self.t50 < 0 and len(self.seen) >= self.converging_at/2:
                    self.t50=self.cycle
                    self.t50_gene_list=self.gene_list.copy()
                if self.t50 >= 0 and self.AB_convergence < 0:
                    if self.converged_AtoB >= self.genesA and self.converged_BtoA >= self.genesB:
                        self.AB_convergence=self.cycle
            if show_output:
                # show the progress bar with level_of_convergence completion
                self._print_progress(self.level_of_convergence*100)

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
        for i in range(self.genesA+self.genesB):
            conv+=i
        return conv-1

    def _pick_cut(self, type):
        """
        return a tuple of indices that represent the breakpoints of the next inversion

        type specifies if the cut is for 'inversion' or 'translocation'
        """

        if self.cycle < len(self.inversion_cuts):
            return self.inversion_cuts[self.cycle]

        valid=False
        while not valid:
            cut=(random.randint(1, len(self.gene_list)-1), random.randint(1, len(self.gene_list)-1))
            valid=self._validate_cut(cut, type)
        return cut

    def _validate_cut(self, cut, type):
        """
        check whether the two breakpoints chosen should be rejected or not
        type specifies if the cut is for 'inversion' or 'translocation'

        no rules are imposed currently, returns True
        """
        
        return True
    
        import random as r
        distance=abs(cut[1]-cut[0])
        if distance == 0:
            return False
        cutoff=1/distance
        rand=r.random()
        return rand <= cutoff

    def _transpose_genes(self):
        """
        translocate translocations_per_cycle genes within the gene list
        """

        for i in range(self.translocations_per_cycle):
            i0, i1=self.pick_cut()
            gene=self.gene_list.pop(i0)
            self.gene_list.insert(i1, gene)
            
    def _shuffle(self):
        """
        invert chromosome and translocate genes

        translocations are not currently done, as they are not implemented with the new _update_cycle process
        """

        self._invert()
        # self._transpose_genes()
        self._update_cycle(self.inversion_cuts[-1])

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
        
    def _update_cycle(self, cut):
        """
        prepare and execute updates of seen and m (only in the window around the breakpoints of cut)
        """

        i0 = cut[0]
        i1 = cut[1]

        # increment cycle
        self.cycle+=1

        # update seen around i0 and i1
        start, end=self._get_window(i0)
        self._update_seen(start, end)
        start, end=self._get_window(i1)
        self._update_seen(start, end)

        # update m around i0 and i1
        self._update_m(i0, i1)
        
    def _update_seen(self, start, end):
        """
        update the seen list and trace structures
        """

        # iterate over all genes between start and end
        for gene_index in range(start, end):
            # iterate over all following genes to compare to the current one
            for comparing_gene_index in range(gene_index+1, min(end+1, gene_index+self.window_size+1)):
                # group the two genes in a tuple
                pair=tuple(sorted([self.gene_list[gene_index], self.gene_list[comparing_gene_index]]))
                # update data structures if the pairing is new
                if pair not in self.seen:
                    # udpate seen with the current cycle
                    self.seen[pair]=self.cycle

                    # update trace structure for both genes in the tuple
                    for pair_index in [0,1]:
                        other=1 if pair_index == 0 else 0

                        # update the current gene's trace with the current cycle and increment the last entry by 1 => (cycle, number_of_interactions)
                        self.trace[pair[pair_index]].append((self.cycle, self.trace[pair[pair_index]][-1][1]+1))

                        # update inter-group traces with a new interaction
                        if pair[pair_index].startswith("A") and pair[other].startswith("B"):
                            self.trace_AtoB[pair[pair_index]].append((self.cycle, self.trace_AtoB[pair[pair_index]][-1][1]+1)) #.[-1] += 1
                            # check whether the latest count of cross-group gene interactions for this gene is equal
                            #   to the number of genes in the opposite group, increase counter for converged A-to-B genes
                            # subtract 1 from the needed interactions if this gene is A.1 (telomeres stay intact and
                            #   cannot interact with the other telomere)
                            if self.trace_AtoB[pair[pair_index]][-1][1] == self.genesB-(1 if pair[pair_index].split('.')[1] == '1' else 0):
                                self.converged_AtoB+=1
                        if pair[pair_index].startswith("B") and pair[other].startswith("A"):
                            self.trace_BtoA[pair[pair_index]].append((self.cycle, self.trace_BtoA[pair[pair_index]][-1][1]+1)) #[-1] += 1
                            # same as above, but for B-to-A
                            if self.trace_BtoA[pair[pair_index]][-1][1] == self.genesA-(1 if pair[pair_index].split('.')[1] == str(self.genesB) else 0):
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
