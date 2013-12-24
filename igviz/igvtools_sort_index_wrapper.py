#!/usr/bin/env python

# Import modules
import subprocess, sys, os, optparse
script_path = os.path.dirname(__file__)
myutils_path = os.path.join(script_path, '../myutils')
sys.path.append(myutils_path)
import myos
import doe_reader
import execs_commands 

def run_igvtools_index(in_fn, index_it, sort_fn, logs_dir):
    qname = 'regevlab'
    mem_usage = '5000'
    job_name = myos.basename_no_ext(in_fn)
    bsubcmd = myos.create_bsub_string_no_rm_logs_dir(logs_dir, job_name, qname = qname, mem_usage = mem_usage)
    if index_it:
        if sort_fn is not None:
            sort_cmd = execs_commands.igvtools().sort(in_fn, sort_fn)
            index_cmd = execs_commands.igvtools().index(sort_fn)
            runcmd = sort_cmd+';'+index_cmd
        else:
            index_cmd = execs_commands.igvtools().index(in_fn)
    else:
        runcmd = execs_commands.igvtools().sort(in_fn, sort_fn)
    fullcmd = bsubcmd+'\"'+runcmd+'\"'
    print fullcmd
    #os.system(fullcmd)
    return 0

def main():
    parser = optparse.OptionParser("Wrapper for igvtools index and sort software -x means index, -x -s means sort then index, -s means sort")
    parser.add_option('-i', '--in_fn', action = "store", help = "igvtools index accepts: .sam, .aligned, .vcf, .psl, and .bed")
    parser.add_option('-x', '--index_it', action = "store_true", default = False, help = "type -x or --index_it If want to index the file")
    parser.add_option('-s', '--sort_fn', action = "store", default = None, help = "igvtools sort accepts: .cn, .igv, .sam, .aligned, .vcf, .psl, and .bed")
    parser.add_option('-l', '--logs_dir', action = "store")
    (options, args) = parser.parse_args()
    run_igvtools_index(options.in_fn, options.index_it, options.sort_fn, options.logs_dir)
    return 0

if __name__ == "__main__":
    main()
