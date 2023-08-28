#!/usr/bin/env python3

"""
This is the file to be executed in in order to simulate inversions on a chromosome.
"""

from chromosome import Chrom
import utils
import plot
import time
import os

def main():
    # handle command line arguments
    parser=utils.create_parser()
    namespace=parser.parse_args()
    args=vars(namespace)
    
    # start a timer
    start=time.time()
    
    outdir=args['output_dir']
    if not os.path.exists(outdir):
        raise parser.error("The directory {} does not exist.".format(outdir))
    if not outdir[-1] == '/':
        outdir+='/'

    if args['simulate']:
        Asize=args['asize']
        Bsize=args['bsize']
        if Asize <= 0 or Bsize <= 0:
            raise parser.error("Asize and Bsize have to be specified and be > 0")
        wsize=args['window_size']
        if not wsize in range(1, Asize+Bsize+1):
            raise parser.error("window-size has to be between 1 and the chromosome size (Asize+Bsize)")
        loc=args['level_of_convergence']
        if not loc in range(1, Asize+Bsize+1):
            raise parser.error("level-of-convergence has to be between 0 and 1")
        cycles=args['cycle_number']
        tpc=args['translocations_per_cycle']
        outname=args['filename']

        chrom=Chrom(Asize, Bsize, level_of_convergence=loc, window_size=wsize, translocations_per_cycle=tpc)
        chrom.run(n=cycles, show_output=True, trace=False)
        utils.save(chrom, outdir, outname)
    
    if args['plot']:
        source=args['source']
        gif=args['gif']
        plot.plot_chrom(source, outdir, gif)
        
    end=time.time()
    elapsed=end-start
    elapsed_string="{minutes:02d}:{seconds:02d}".format(minutes=int(elapsed//60), seconds=int(elapsed%60))
    print("\nelapsed time: "+elapsed_string)
    
if __name__ == "__main__":
    main()
