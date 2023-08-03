#!/usr/bin/env python3

from chromosome import Chrom
import utils
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
    converge=args['converge']
    cycles=args['cycle_number']
    if not converge and cycles <= 0:
        raise parser.error("please specify a cycle number >0 or use --converge")

    tpc=args['translocations_per_cycle']

    do_plots=args['plot_curves']
    
    #do_average_t50=args['average_t50']
    #if do_average_t50:
    #    average=utils.calculate_average_t50(outdir, output_name)
    #    print("average t50: {average:.2f} cycles".format(average=average))
    #    return
    
    # start a timer
    start=time.time()

    print("\ncreating chromosome...")
    chrom=Chrom(0, Asize, Bsize, level_of_convergence=loc, window_size=wsize, until_converged=converge, translocations_per_cycle=tpc)
    print("\nrunning simulation...")
    chrom.simulation_cycle(iterations=cycles)
    if do_plots:
        print("\nplotting results...")
        utils.plot_results(chrom, outdir)

    end=time.time()
    elapsed=end-start
    elapsed_string="{minutes:02d}:{seconds:02d}".format(minutes=int(elapsed//60), seconds=int(elapsed%60))
    print("\nelapsed time: "+elapsed_string)

    utils.log(chrom, outdir, elapsed=elapsed_string)
    
if __name__ == "__main__":
    main()
