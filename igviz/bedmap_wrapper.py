#!/usr/bin/env python

# Import modules
import subprocess, sys, os, optparse
script_path = os.path.dirname(__file__)
myutils_path = os.path.join(script_path, '../myutils')
sys.path.append(myutils_path)
import myos
import doe_reader
import execs_commands 

def run_bedmap(in_reference_fn, doe_csv, header_name_of_exp_id, in_map_dir, extension, bedmap_options, out_dir, logs_dir):
    qname = 'regevlab'
    mem_usage = '5000'
    bed_map_fn_dict = doe_reader.create_experiment_fns(doe_csv, header_name_of_exp_id, in_map_dir, extension)
    for exp_name, bed_map_fn in bed_map_fn_dict.iteritems():
        bsubcmd = myos.create_bsub_string_no_rm_logs_dir(logs_dir, myos.basename_no_ext(bed_map_fn), qname = qname, mem_usage = mem_usage)
        out_fn = os.path.join(out_dir, myos.basename_no_ext(bed_map_fn)+'_mapped2_'+myos.basename_no_ext(in_reference_fn)+'.bedmap')
        runcmd = execs_commands.bedops().bedmap(bedmap_options, in_reference_fn, bed_map_fn, out_fn)
        fullcmd = bsubcmd+'\"'+runcmd+'\"'
        print fullcmd
        os.system(fullcmd)
    return 0

def main():
    parser = optparse.OptionParser("Wrapper for phasedBam2bed nimrod script that divides a bam into two pat/mat beds")
    parser.add_option('-r', '--in_reference_fn', action = "store")
    parser.add_option('-d', '--doe_csv_fn', action = "store")
    parser.add_option('--header_name_of_exp_id', action = "store")
    parser.add_option('-i', '--in_map_dir', action = "store")
    parser.add_option('-e', '--extension_after_exp_id', action = "store")
    parser.add_option('--bedmap_options', help = "--bedmap_options='--faster --count --echo'", action = "store")
    parser.add_option('-o', '--out_dir')
    parser.add_option('-l', '--logs_dir', action = "store")
    (options, args) = parser.parse_args()
    run_bedmap(options.in_reference_fn, options.doe_csv_fn, options.header_name_of_exp_id, options.in_map_dir, options.extension_after_exp_id, options.bedmap_options, options.out_dir, options.logs_dir)
    return 0

if __name__ == "__main__":
    main()
