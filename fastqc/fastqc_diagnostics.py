#!/usr/bin/env python

# Import modules
import subprocess, sys, os, optparse
script_path = os.path.dirname(__file__)
myutils_path = os.path.join(script_path, '../myutils')
sys.path.append(myutils_path)
import myos
import doe_reader
SCRIPTS = '/seq/epiprod/de/scripts/nimrod/'
dependencies_list = ['.fastqc-0.10.1']

def run_fastqc(doe_csv_fn, in_dir, extension_name, logs_dir):
    qname = 'regevlab'
    mem_usage = '5000'
    dict_fq_fns = doe_reader.create_experiment_fns(doe_csv_fn, 'name', in_dir, extension_name)
    out_dir = in_dir
    dep_cmd = myos.load_dependencies_cmd(dependencies_list)
    dep_cmd = dep_cmd+'reuse .fastqc-0.10.1;'
    for exp_name, in_fn in dict_fq_fns.iteritems():
        bsubcmd = myos.create_bsub_string_no_rm_logs_dir(logs_dir, exp_name, qname = qname, mem_usage = mem_usage)
        runcmd = myos.fastqc(in_fn, out_dir)
        fullcmd = bsubcmd+' \"'+runcmd+'\"'
        print fullcmd
        #os.system(fullcmd)
    return 0

def main():
    parser = optparse.OptionParser()
    parser.add_option('-d', '--doe_csv_fn', action = "store") 
    parser.add_option('-i', '--in_dir', action = "store") 
    parser.add_option('-l', '--logs_dir', action = "store", help = "a log dir or n for default out_dir+'trimmed/'")
    parser.add_option('-e', '--extension_name', action = "store", help = "i.e., .fastq, .gz, .fq, _trimmed.fastqc, etc. (add the dot if it's a .extn, or add the _ if it's a _ extn) it will run fastqc for each file with such extension in the in_dir; for a single file just do file_name, i.e, MF1i5_P8_trimmed.fastq")
    myargs, noidea = parser.parse_args()
    run_fastqc(myargs.doe_csv_fn, myargs.in_dir, myargs.extension_name, myargs.logs_dir)
    return 0

if __name__ == "__main__":
    main()
