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
            exp_name = row[0].split(',')[3]
            adapter_seq = row[0].split(',')[5]
            ens[exp_name] = [adapter_seq, os.path.join(in_dir, exp_name+'.'+extension_name)]
    return ens

def run_fastx_fastq(logs_dir, trim_options, ADf_options, Qf_options, Af_options, in_dir, out_dir, doe_csv, extension_name):
    wlist_out = []
    qname = 'hour'
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
        out_Cf_fn = os.path.join(out_dir, exp_name+'.cf.fastq')
        runcmd_Cf = casava_quality_filter(in_fn, out_Cf_fn)
        out_ADf_fn = os.path.join(out_dir, exp_name+'.cf.adf.fastq')
        runcmd_ADf = fastx_adaptor_filter(adaptor_seq, ADf_options, out_Cf_fn, out_ADf_fn)
        rm_Cf_fn = 'rm %s;' %(out_Cf_fn)
        out_Qf_fn = os.path.join(out_dir, exp_name+'.cf.adf.qf.fastq')
        runcmd_Qf = fastq_quality_filter(Qf_options, out_ADf_fn, out_Qf_fn)
        rm_ADf_fn = 'rm %s;' %(out_ADf_fn)
        out_Af_fn = os.path.join(out_dir, exp_name+'.cf.adf.qf.af.fastq')
        runcmd_Af = fastx_artifacts_filter(Af_options, out_Qf_fn, out_Af_fn)
        rm_Qf = 'rm %s;' %(out_Qf_fn)
        out_trim_fn = os.path.join(out_dir, exp_name+'.cf.adf.qf.af.trim.fastq')
        runcmd_trim = fastx_trimmer(trim_options, out_Af_fn, out_trim_fn)
        rm_Af = 'rm %s;' %(out_Af_fn)
        fullcmd = bsubcmd+' \"'+''.join([dep_cmd,runcmd_Cf,runcmd_ADf,rm_Cf_fn,runcmd_Qf,rm_ADf_fn,runcmd_Af,rm_Qf,runcmd_trim])+'\"'
        #fullcmd = bsubcmd+' \"'+''.join([dep_cmd,runcmd_Cf,runcmd_ADf,rm_Cf_fn,runcmd_Qf,rm_ADf_fn,runcmd_Af,rm_Qf,runcmd_trim,rm_Af])+'\"'
        print fullcmd
        #os.system(fullcmd)
    return 0

def run_fastqc(logs_dir, in_dir):
    qname = 'hour'
    mem_usage = '5000'
    job_name = 'qc'
    out_log_fn = os.path.join(logs_dir, job_name+'.out')
    err_log_fn = os.path.join(logs_dir, job_name+'.err')
    bsubcmd = create_bsub_string_no_remove(logs_dir, job_name, out_log_fn, err_log_fn, qname, mem_usage)
    in_fns = glob.glob('%s/*.cf' %(in_dir))
    in_fns = [os.path.basename(fn) for fn in in_fns]
    in_fns = ' '.join(in_fns)
    out_dir = in_dir
    #dep_cmd = load_dependencies_cmd(dependencies_list)
    #dep_cmd = dep_cmd+'reuse .fastqc-0.10.1;'
    runcmd = fastqc(in_fns, out_dir)
    fullcmd = bsubcmd+' \"'+runcmd+'\"'
    print fullcmd
    os.system(fullcmd)
    return 0


def main():
    parser = optparse.OptionParser()
    parser.add_option('--trim_options', action = "store", help = "options as a string i.e., --options='-Q 33 -f 16 -v'")
    parser.add_option('--Qf_options', action = "store", help = "options as a string i.e., --options='-Q 33 -q 20 -p 90 -v'")
    parser.add_option('--Af_options', action = "store", help = "options as a string i.e., --options='-Q 33 -v'")
    parser.add_option('--ADf_options', action = "store", help = "options as a string i.e., --options='-Q 33 -C', -a ATTCGC added at the beginning")
    parser.add_option('-q', '--fastqc', action = "store", help = "y/n to also run fastqc")
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
        out_dir = os.path.join(in_dir,'trimmed')
    if logs_dir == 'n':
        logs_dir = os.path.join(out_dir, 'logs')
    run_fastx_fastq(logs_dir, myargs.trim_options, myargs.ADf_options, myargs.Qf_options, myargs.Af_options, myargs.in_dir, out_dir, myargs.doe_csv, extension_name)
    if myargs.fastqc == 'y':
        run_fastqc(logs_dir, myargs.in_dir)
    return 0


if __name__ == "__main__":
    main()
