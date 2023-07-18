from chromosome import Chrom
import utils
import time

def main():
    # handle command line arguments
    parser=utils.get_parser
    namespace=parser.parse_args()
    args=vars(namespace)
    outdir=args['output_dir']
    if not os.path.exists(outdir):
        raise parser.error("The directory {} does not exist.".format(outdir))
    if not outdir[-1] == '/':
        outdir+='/'
    
    # start a timer
    start=time.time()

    print("\ncreating chromosome...")
    chrom=Chrom(0, args['Asize'], args['Bsize'], level_of_convergence=args['level_of_convergence'])
    print("\nrunning simulation...")
    chrom.simulation_cycle(until_converged=args['converge'])
    print("\nplotting results...")
    utils.plot_results(chrom, outdir)

    end=time.time()
    elapsed=end-start
    elapsed_string="{minutes:02d}:{seconds:02d}".format(minutes=int(elapsed//60), seconds=int(elapsed%60))
    print("\nelapsed time: "+elapsed_string)
    
    utils.log(chrom, outdir, elapsed=elapsed_string)
    
if __name__ == "__main__":
    main()
