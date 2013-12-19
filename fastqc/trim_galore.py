#!/usr/bin/env python

# Import modules
import subprocess, sys, os, glob, optparse, csv
sys.path.append('/seq/epigenome01/allelix/rusty/auto_homer/myos/')
from myos import *
#SCRIPTS = '${de}scripts/nimrod/'
SCRIPTS = '/seq/epiprod/de/scripts/nimrod/'
dependencies_list = ['.fastx-toolkit-0.0.13', '.fastqc-0.10.1']

def read_experiment_names(doe_csv, in_dir, extension_name):
    ens = {}
    with open(doe_csv,'rb') as csvfile:
        reader = csv.reader(csvfile, delimiter=' ', quotechar='|')
        for row in reader:
            row_list = row[0].split(',')
            exp_name = row_list[3]
            if exp_name == 'name':
                continue
            adapter_seq = row_list[7]  # here we take the adaptorStarPlusIndex
            seq_length = str(int(row_list[4])-14-6)
            ens[exp_name] = [adapter_seq, os.path.join(in_dir, exp_name+'.'+extension_name), seq_length]
    print ens
    return ens

def run_trim_galore(logs_dir, in_dir, out_dir, doe_csv, trim_galore_options, extension_name):
    wlist_out = []
    qname = 'regevlab'
    mem_usage = '5000'
    exp_names = read_experiment_names(doe_csv, in_dir, extension_name) # experiment names
    #in_fns = glob.glob('%s/*.cf' %(in_dir))
    dep_cmd = load_dependencies_cmd(dependencies_list)
    index = 0
    remove_all_files_given_dir(out_dir)
    check_if_directory_exists_create_it(out_dir)
    for exp_name in exp_names:
        index+=1
        adaptor_seq = exp_names[exp_name][0]
        in_fn = exp_names[exp_name][1]
        job_name = exp_name
        out_log_fn = os.path.join(logs_dir, job_name+'.out')
        err_log_fn = os.path.join(logs_dir, job_name+'.err')
        bsubcmd = create_bsub_string_no_remove(logs_dir, job_name, out_log_fn, err_log_fn, qname, mem_usage)
        runcmd_tgf = trim_galore_filter(adaptor_seq, trim_galore_options+' --length %s' %(exp_names[exp_name][2]), in_fn, out_dir)
        fullcmd = bsubcmd+' \"'+''.join([dep_cmd,runcmd_tgf])+'\"'
        print fullcmd
        os.system(fullcmd)
    return 0

def main():
    parser = optparse.OptionParser()
    parser.add_option('--trim_galore_options', action = "store", help = "options as a string i.e., --trim_galore_options='--clip_R1 10 --clip_R2 10 --gzip --phred33 -q 15 -e 0.05 -s 5 --length 48'")
    #parser.add_option('-c', '--casava_filter', action = "store", help = "y/n to also run casava_filter")
    #parser.add_option('-q', '--quality_filter', action = "store", help = "y/n to also filter by quality")
    #parser.add_option('-a', '--artifacts_filter', action = "store", help = "y/n to also filter artifacts")
    #parser.add_option('-d', '--adaptor_filter', action = "store", help = "y/n to also filter reads that contain the adaptor sequence")
    parser.add_option('-i', '--in_dir', action = "store") 
    parser.add_option('-o', '--out_dir', action = "store", help = "an out dir or n for default in_dir+'trimmed/'")
    parser.add_option('-l', '--logs_dir', action = "store", help = "a log dir or n for default out_dir+'trimmed/'")
    parser.add_option('-d', '--doe_csv', action = "store", help = "fullname to csv files containing experiment names and adaptor seqs")
    myargs, noidea = parser.parse_args()
    in_dir = myargs.in_dir
    extension_name = 'fastq'
    out_dir = myargs.out_dir
    logs_dir = myargs.logs_dir
    if out_dir == 'n':
        out_dir = os.path.join(in_dir,'trim_galore')
    if logs_dir == 'n':
        logs_dir = os.path.join(out_dir, 'trim_galore_logs')
    run_trim_galore(logs_dir, in_dir, out_dir, myargs.doe_csv, myargs.trim_galore_options, extension_name)
    return 0

if __name__ == "__main__":
    main()
