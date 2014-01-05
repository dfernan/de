#!/usr/bin/env python

# Import modules
import subprocess, sys, os, optparse
script_path = os.path.dirname(__file__)
myutils_path = os.path.join(script_path, '../myutils')
sys.path.append(myutils_path)
import myos
import doe_reader
import execs_commands 
import general
from itertools import izip

def average(tos, map_lst):
    """
    given
        tos: a sequence of N|user\tname\tAGE\n
        map_lst: a list with positions corresponding to those in tos, and values
                 indicating which group each tos element will be averaged with.
    return the groups of averages as a list of user\tname\tAVG\n
    http://stackoverflow.com/questions/20893768/function-to-compute-averages-and-retrieve-information-over-defined-indexes-in-a/20895444?noredirect=1#20895444
    """
    # get the leading nums
    nums = [s.partition('|')[0] for s in tos]
    # group them into lists that will be averaged together (based on the map)
    avg_groups = [[] for i in set(map_lst)]
    for i,n in zip(map_lst, nums):
        avg_groups[i].append(float(n))
    # generate the averages
    def fmt(tup):
        tmp = tos[0].partition('|')[2]
        tmp = tmp.strip().split('\t')[0:3]
        mid = '\t'.join(tmp)
        avg = round(sum(tup)/float(len(tup)),3)
        return "{0}\t{1}\n".format(mid, avg)
    return [fmt(l) for l in avg_groups]

def close_lst_of_fhs(lst_of_fhs):
    for fh in lst_of_fhs:
        fh.close()
    return 0

def run_bedmap(doe_csv, header_name_of_exp_id, header_names_of_factors, in_map_dir, parent, extension, out_dir):
    qname = 'regevlab'
    mem_usage = '5000'
    bed_map_fn_dict = doe_reader.create_experiment_fns(doe_csv, header_name_of_exp_id, in_map_dir, '.'+parent+extension)
    factors_dict = doe_reader.read_experiment_fields(doe_csv, header_name_of_exp_id, header_names_of_factors.split(','))
    factors_name = [name for name in factors_dict]
    factors_set = set([':'.join(factors_dict[name]) for name in factors_dict])
    factors_name_corresponding_set = [0 for i in range(0,len(factors_name))]
    i = 0
    for name, factors in factors_dict.iteritems():
        factors = ':'.join(factors)
        index_of_factors = general.index_in_unique_list(factors_set, factors)
        factors_name_corresponding_set[i] = index_of_factors
        i+=1
    files = [open(bed_map_fn_dict[name]) for name in bed_map_fn_dict]
    out_files = [open(os.path.join(out_dir, factor+'.%s.bedgraph' %(parent)), 'w') for factor in factors_set]
    for lines in izip(*files):
        lines_to_print = average(lines, factors_name_corresponding_set)
        f_index = 0
        for line_to_print in lines_to_print:
            out_files[f_index].write(line_to_print)
            f_index += 1
    close_lst_of_fhs(files)
    close_lst_of_fhs(out_files)
    return 0

def main():
    parser = optparse.OptionParser("Wrapper for phasedBam2bed nimrod script that divides a bam into two pat/mat beds")
    parser.add_option('-d', '--doe_csv_fn', action = "store")
    parser.add_option('--header_name_of_exp_id', action = "store")
    parser.add_option('--header_names_of_factors', action = "store", help = "separated by comma, i.e., --header_names_of_factors='age,cross' factors to combine for the average")
    parser.add_option('-i', '--in_map_dir')
    parser.add_option('-p', '--parent', action = "store", help = "parent, m or p")
    parser.add_option('-e', '--extension_after_exp_id', action = "store", help = "extension after exp_id.parent")
    parser.add_option('-o', '--out_dir')
    (options, args) = parser.parse_args()
    run_bedmap(options.doe_csv_fn, options.header_name_of_exp_id, options.header_names_of_factors, options.in_map_dir, options.parent, options.extension_after_exp_id, options.out_dir)
    return 0

if __name__ == "__main__":
    main()
