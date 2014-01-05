#!/usr/bin/env python

# Import modules
import subprocess, sys, os, optparse
script_path = os.path.dirname(__file__)
myutils_path = os.path.join(script_path, '../myutils')
sys.path.append(myutils_path)
import myos
import general
import doe_reader
import execs_commands 

def run_cp_chris2rawdata(logs_dir, out_dir, doe_csv_fn, header_name_of_in_fn):
    qname = 'regevlab'
    mem_usage = '5000'
    dict_fq_fns = doe_reader.read_experiment_field(doe_csv_fn, 'name', header_name_of_in_fn)
    myos.remove_all_files_given_dir(out_dir)
    myos.check_if_directory_exists_create_it(out_dir)
    for exp_name, in_fn in dict_fq_fns.iteritems():
        bsubcmd = myos.create_bsub_string_no_rm_logs_dir(logs_dir, exp_name, qname = qname, mem_usage = mem_usage)
        if os.path.splitext(in_fn)[1] == '.gz':
            cp_cmd = 'cp %s %s; gunzip %s' %(in_fn, os.path.join(out_dir, exp_name+'.fq.gz'), os.path.join(out_dir, exp_name+'.fq.gz'))
        else:
            cp_cmd = 'cp %s %s' %(in_fn, os.path.join(out_dir, exp_name+'.fq'))
        fullcmd = bsubcmd+'\"'+cp_cmd+'\"'
        print fullcmd
        os.system(fullcmd)
    return 0

def main():
    parser = optparse.OptionParser("Input: fastq files in DoE csv file, Output: trimmed fastq files with trim_galore")
    parser.add_option('-o', '--out_dir', action = "store")
    parser.add_option('-l', '--logs_dir', action = "store") 
    parser.add_option('-d', '--doe_csv', action = "store", help = "fullname to csv files containing experiment names and adaptor seqs")
    parser.add_option('-n', '--header_name_of_in_fn', default = None, action = "store")
    (options, args) = parser.parse_args()
    run_cp_chris2rawdata(options.logs_dir, options.out_dir, options.doe_csv, options.header_name_of_in_fn)
    return 0

if __name__ == "__main__":
    main()
