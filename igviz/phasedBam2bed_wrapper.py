#!/usr/bin/env python

# Import modules
import subprocess, sys, os, optparse
script_path = os.path.dirname(__file__)
myutils_path = os.path.join(script_path, '../myutils')
sys.path.append(myutils_path)
import myos
import doe_reader
import execs_commands 

def run_phasedBam2bed(doe_csv, in_dir, out_dir, header_name_of_exp_id, extension, logs_dir):
    qname = 'regevlab'
    mem_usage = '5000'
    if extension is not None:
        bam_fn_dict = doe_reader.create_experiment_fns(doe_csv, header_name_of_exp_id, in_dir, extension)
    else:
        bam_fn_dict = doe_reader.create_experiment_fns(doe_csv, header_name_of_exp_id, in_dir, '.bam')
    for exp_name, bam_in_fn in bam_fn_dict.iteritems():
        bsubcmd = myos.create_bsub_string_no_rm_logs_dir(logs_dir, exp_name, qname = qname, mem_usage = mem_usage)
        bed_p_fn = os.path.join(out_dir, exp_name+'.p.bed')
        bed_m_fn = os.path.join(out_dir, exp_name+'.m.bed')
        runcmd = execs_commands.phasedBam2bed(bam_in_fn, bed_p_fn, bed_m_fn)
        fullcmd = bsubcmd+'\"'+runcmd+'\"'
        print fullcmd
        os.system(fullcmd)
    return 0

def main():
    parser = optparse.OptionParser("Wrapper for phasedBam2bed nimrod script that divides a bam into two pat/mat beds")
    parser.add_option('-d', '--doe_csv_fn', action = "store")
    parser.add_option('-i', '--in_dir_bam_files', action = "store")
    parser.add_option('-o', '--out_dir', action = "store", help = 'outpout directory where to output the .p.bed and .m.bed files')
    parser.add_option('--header_name_of_exp_id', action = "store")
    parser.add_option('-e', '--extension_after_exp_id', default = None, action = "store")
    parser.add_option('-l', '--logs_dir', action = "store")
    (options, args) = parser.parse_args()
    run_phasedBam2bed(options.doe_csv_fn, options.in_dir_bam_files, options.out_dir, options.header_name_of_exp_id, options.extension_after_exp_id, options.logs_dir)
    return 0

if __name__ == "__main__":
    main()
