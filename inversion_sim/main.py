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

    outdir=args['output_dir']
    if not os.path.exists(outdir):
        raise parser.error("The directory {} does not exist.".format(outdir))
    if not outdir[-1] == '/':
        outdir+='/'

    logmeta=args['metalog']
    if logmeta:
        utils.metalog(outdir)
    avt50=args['plot_average_t50']
    if avt50:
        utils.plot_average_t50s(outdir)
    if avt50 or logmeta:
        return
        
    Asize=args['Asize']
    Bsize=args['Bsize']
    if not Asize or not Bsize:
        raise parser.error("Asize and Bsize have to be specified and be >0")
    
    wsize=args['window_size']
    if not wsize in range(1, Asize+Bsize+1):
        raise parser.error("window-size has to be between 1 and the chromosome size (Asize+Bsize)")

    loc=args['level_of_convergence']
    cycles=args['cycle_number']
    
    tpc=args['translocations_per_cycle']

    do_plots=args['plot_curves']
        
    # start a timer
    start=time.time()

    print("\ncreating chromosome...")
    chrom=Chrom(Asize, Bsize, level_of_convergence=loc, window_size=wsize, translocations_per_cycle=tpc)
    print("\nrunning simulation...")
    chrom.run(n=cycles)
    if do_plots:
        print("\nplotting results...")
        plot.plot_results(chrom, outdir)
        
    #utils.log(chrom, outdir, elapsed=elapsed_string)
    utils.save(chrom)
        
    end=time.time()
    elapsed=end-start
    elapsed_string="{minutes:02d}:{seconds:02d}".format(minutes=int(elapsed//60), seconds=int(elapsed%60))
    print("\nelapsed time: "+elapsed_string)
    
if __name__ == "__main__":
    main()
