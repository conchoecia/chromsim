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

    Asize=args['Asize']
    Bsize=args['Bsize']
    if Asize <= 0 or Bsize <= 0:
        raise parser.error("Asize and Bsize have to be > 0")
    
    wsize=args['window_size']
    if not wsize in range(1, Asize+Bsize+1):
        raise parser.error("window-size has to be between 1 and the chromosome size (Asize+Bsize)")

    loc=args['level_of_convergence']
    converge=args['converge']
    
    output_name='inversion_sim_A{Asize}-B{Bsize}_l{loc}_w{wsize}'.format(Asize=Asize, Bsize=Bsize, loc=loc, wsize=wsize)+('_c' if converge else '')
    
    # start a timer
    start=time.time()

    print("\ncreating chromosome...")
    chrom=Chrom(0, Asize, Bsize, level_of_convergence=loc, window_size=wsize)
    print("\nrunning simulation...")
    chrom.simulation_cycle(until_converged=converge)
    print("\nplotting results...")
    utils.plot_results(chrom, outdir, output_name)

    end=time.time()
    elapsed=end-start
    elapsed_string="{minutes:02d}:{seconds:02d}".format(minutes=int(elapsed//60), seconds=int(elapsed%60))
    print("\nelapsed time: "+elapsed_string)

    log_file_name='inversion_sim'
    utils.log(chrom, outdir, output_name, elapsed=elapsed_string)
    
if __name__ == "__main__":
    main()
