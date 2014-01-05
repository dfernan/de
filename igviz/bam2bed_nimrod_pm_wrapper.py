#!/usr/bin/env python

# Import modules
import subprocess, sys, os, optparse
script_path = os.path.dirname(__file__)
myutils_path = os.path.join(script_path, '../myutils')
sys.path.append(myutils_path)
import myos
import doe_reader
import execs_commands 

def run_vcf2bed(vcf_in_fn, vcf_out_fn, logs_dir):
    qname = 'regevlab'
    mem_usage = '5000'
    job_name = myos.basename_no_ext(vcf_in_fn)
    bsubcmd = myos.create_bsub_string_no_rm_logs_dir(logs_dir, job_name, qname = qname, mem_usage = mem_usage)
    runcmd = execs_commands.bedops().vcf2bed(vcf_in_fn, vcf_out_fn)
    fullcmd = bsubcmd+'\"'+runcmd+'\"'
    print fullcmd
    #os.system(fullcmd)
    return 0

def main():
    parser = optparse.OptionParser("Wrapper for BedOps vcf2bed script")
    parser.add_option('-i', '--vcf_in_fn', action = "store")
    parser.add_option('-o', '--vcf_out_fn', action = "store")
    parser.add_option('-l', '--logs_dir', action = "store")
    (options, args) = parser.parse_args()
    run_vcf2bed(options.vcf_in_fn, options.vcf_out_fn, options.logs_dir)
    return 0

if __name__ == "__main__":
    main()
