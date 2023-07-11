import os
import matplotlib.pyplot as plt

def convert_dt(dt):
    mins, secs=dt.split(':')
    return 60*int(mins)+int(secs)

def main():
    remote_path='/scratch/molevo/bluehmel/'
    use_dir=(remote_path if os.path.exists(remote_path) else '')+'log/'
    input_file='log.csv'
    output_file='runtime.png'

    if not os.path.exists(use_dir+input_file):
        print(use_dir+input_file)
        print('log file does not exist')
        return

    data=None
    
    with open(use_dir+input_file) as f:
        lines=[line.rstrip().split(';') for line in f]
        header=lines[0]
        pos_Asize=header.index('|A|')
        pos_Bsize=header.index('|B|')
        pos_dt=header.index('Delta_t')
        data=[(int(line[pos_Asize])+int(line[pos_Bsize]), convert_dt(line[pos_dt])) for line in lines[1:]]

    plt.plot([x[0] for x in data], [x[1]/60 for x in data], 'o')
    plt.xlabel("chromosome size (|A|+|B|)")
    plt.ylabel("runtime (minutes)")
    plt.savefig(use_dir+output_file)

if __name__ == '__main__':
    main()
