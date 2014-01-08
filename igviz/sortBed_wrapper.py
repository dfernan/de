#!/usr/bin/env python

# Import modules
import subprocess, sys, os, optparse
script_path = os.path.dirname(__file__)
myutils_path = os.path.join(script_path, '../myutils')
sys.path.append(myutils_path)
import myos
import doe_reader
import execs_commands 

def run_sortBed(doe_csv, in_dir, header_name_of_exp_id, extension, logs_dir, remove, execute):
    qname = 'regevlab'
    mem_usage = '20000'
    bed_fn_dict = doe_reader.create_experiment_fns(doe_csv, header_name_of_exp_id, in_dir, extension)
    for exp_name, bed_in_fn in bed_fn_dict.iteritems():
        bsubcmd = myos.create_bsub_string_no_rm_logs_dir(logs_dir, exp_name, qname = qname, mem_usage = mem_usage)
        sort_cmd = execs_commands.bedops().sortbed(bed_in_fn, os.path.splitext(bed_in_fn)[0]+'.sorted.bed', '')
        if remove:
            rm_m_cmd = 'rm %s' %(bed_m_fn)
            runcmd = sort_cmd+';'+rm_cmd
        else:
            runcmd = sort_cmd
        fullcmd = bsubcmd+'\"'+runcmd+'\"'
        print fullcmd
        if execute:
            os.system(fullcmd)
    return 0

def main():
    parser = optparse.OptionParser("Wrapper for sortBed")
    parser.add_option('-d', '--doe_csv_fn', action = "store")
    parser.add_option('-i', '--in_dir_bam_files', action = "store")
    parser.add_option('--header_name_of_exp_id', action = "store")
    parser.add_option('-e', '--extension_after_exp_id', default = None, action = "store")
    parser.add_option('-l', '--logs_dir', action = "store")
    parser.add_option('-r', '--rm', action = "store_true", default=False)
    parser.add_option('-x', '--execute', action = "store_true", default=False)
    (options, args) = parser.parse_args()
    run_sortBed(options.doe_csv_fn, options.in_dir_bam_files, options.header_name_of_exp_id, options.extension_after_exp_id, options.logs_dir, options.rm, options.execute)
    return 0

if __name__ == "__main__":
    main()
