# This module is used by:
# homer_pipeline_run.py, ...
import glob
import os
import pysam
import sys
import random
import re
import subprocess

# Define constants
#HEP = "/seq/epigenome01/Rusty/Homer/homer/bin/" # Homer executables path 
#PEP = "/seq/software/picard/current/bin/" # Picard executables path
#ITEN = "/xchip/igv/tools/./igvtools" # IGV tools executable fullname
TSEP='/seq/epiprod/de/software/tuxsim/src/'


# Define functions


def exit(self, status=0, message=None):
  if message:
    self._print_message(message, _sys.stderr)
  sys.exit(status)


def check_if_directory_exists_create_it(dir):
  if not os.path.exists(dir):
    try:
      os.makedirs(dir)
    except:
      print('Could not create the dir %s' %(dir))
  return 0


def erase_file_if_exists(fullname):
  try:
    os.remove(fullname)
  except OSError:
    pass
  return 0


def remove_all_files_given_dir(dir):
  ''' remove all the files in a given dir, but  keep the directory '''
  dir_plus_star = os.path.join(dir, '*')
  try:
    files = glob.glob(dir_plus_star)
    for f in files:
      os.remove(f)
  except:
    pass
  return 0


def create_bsub_string_rm_logs_dir(logs_dir, job_name, out_log_fn='', err_log_fn='', qname = 'week', mem_usage = 'default'):
    if out_log_fn == '':
        out_log_fn = os.path.join(logs_dir, job_name+'.out')
    if err_log_fn == '':
        err_log_fn = os.path.join(logs_dir, job_name+'.err')
    if os.popen('echo $USER').read().strip() == 'dfernand':
        if mem_usage == 'default':
            s = 'bsub -P epigenome -q %s -J %s -o %s -e %s' %(qname, job_name, out_log_fn, err_log_fn)
        else:
            s = 'bsub -P epigenome -R "rusage[mem=%s]" -q %s -J %s -o %s -e %s' %(mem_usage, qname, job_name, out_log_fn, err_log_fn)
    elif os.popen('echo $USER').read().strip() == 'fernandez':
        s = 'sbatch -J %s -o %s -e %s -p general' %(job_name, out_log_fn, err_log_fn)
    check_if_directory_exists_create_it(logs_dir)
    remove_all_files_given_dir(logs_dir)
    erase_file_if_exists(out_log_fn)
    erase_file_if_exists(err_log_fn)
    return s


def create_bsub_string_no_rm_logs_dir(logs_dir, job_name, out_log_fn='', err_log_fn='', qname = 'week', mem_usage = 'default'):
    if out_log_fn == '':
        out_log_fn = os.path.join(logs_dir, job_name+'.out')
    if err_log_fn == '':
        err_log_fn = os.path.join(logs_dir, job_name+'.err')
    if mem_usage == 'default':
        s = 'bsub -P epigenome -q %s -J %s -o %s -e %s' %(qname, job_name, out_log_fn, err_log_fn)
    else:
        s = 'bsub -P epigenome -R "rusage[mem=%s]" -q %s -J %s -o %s -e %s' %(mem_usage, qname, job_name, out_log_fn, err_log_fn)
    check_if_directory_exists_create_it(logs_dir)
    erase_file_if_exists(out_log_fn)
    erase_file_if_exists(err_log_fn)
    return s


def create_bsub_string_no_remove_matlab(logs_dir, job_name, out_log_fn, err_log_fn, qname = 'week', mem_usage = 'default'):
  if mem_usage == 'default':
    s = 'bsub -P epigenome -R "rusage[matlab=1:duration=1]" -q %s -J %s -o %s -e %s' %(qname, job_name, out_log_fn, err_log_fn)
  else:
    s = 'bsub -P epigenome -R "rusage[matlab=1:duration=1]" -R "rusage[mem=%s]" -q %s -J %s -o %s -e %s' %(mem_usage, qname, job_name, out_log_fn, err_log_fn)
  check_if_directory_exists_create_it(logs_dir)
  erase_file_if_exists(out_log_fn)
  erase_file_if_exists(err_log_fn)
  return s


def create_wstring_from_wlist(wlist):
  by = 'job_id' # or by = 'job_name'
  wstring = '"done('
  for el in wlist:
    if by == 'job_name':
      wstring += '\'%s\') && done(' %(el)
    elif by == 'job_id':
      wstring += '%s) && done(' %(el)
  wstring = wstring[0:-9]+'"'
  return wstring

    
def create_bsub_string_with_dependencies(wstring, logs_dir, job_name, out_log_fn, err_log_fn, qname='week', mem_usage='default'):
  if mem_usage == 'default':
    s = 'bsub -P epigenome -q %s -J %s -o %s -e %s -w %s' %(qname, job_name, out_log_fn, err_log_fn, wstring)
  else:
    s = 'bsub -P epigenome -R rusage[mem=%s] -q %s -J %s -o %s -e %s -w %s' %(mem_usage, qname, job_name, out_log_fn, err_log_fn, wstring)
  check_if_directory_exists_create_it(logs_dir)
  remove_all_files_given_dir(logs_dir)
  erase_file_if_exists(out_log_fn)
  erase_file_if_exists(err_log_fn)
  return s


def create_bsub_string_with_dependencies_matlab(wstring, logs_dir, job_name, out_log_fn, err_log_fn, qname='week', mem_usage='default'):
  if mem_usage == 'default':
    s = 'bsub -P epigenome -R "rusage[matlab=1:duration=1]" -q %s -J %s -o %s -e %s -w %s' %(qname, job_name, out_log_fn, err_log_fn, wstring)
  else:
    s = 'bsub -P epigenome -R "rusage[matlab=1:duration=1]" -R rusage[mem=%s] -q %s -J %s -o %s -e %s -w %s' %(mem_usage, qname, job_name, out_log_fn, err_log_fn, wstring)
  check_if_directory_exists_create_it(logs_dir)
  remove_all_files_given_dir(logs_dir)
  erase_file_if_exists(out_log_fn)
  erase_file_if_exists(err_log_fn)
  return s


def load_dependencies_cmd(dependencies_list):
    dotkit_cmd = ("export DK_ROOT=/broad/software/dotkit;"
                  "source /broad/software/dotkit/ksh/.dk_init;")
    use_cmd = ""
    for dep in dependencies_list:
        use_cmd += 'use %s;' %(dep)
    return dotkit_cmd+use_cmd


def check_if_chr(bam_in_fn):
  ''' this functions checks if a BAM header uses chr reference (true), or no chr reference (false) '''
  s = "samtools view -H %s | grep -o 'SN:c'" %(bam_in_fn)
  a = subprocess.Popen(s,stdout=subprocess.PIPE,shell=True)
  output = a.stdout.read()
  if len(output) == 0:
    return False
  else:
    return True


def create_chromosome_string(chrom_sizes_fn):
  fh = open(chrom_sizes_fn, 'r')
  chrom_list = []
  for line in fh:
    chromosome = line.strip().split('\t')[0]
    a = re.search('(\.|_)',chromosome)
    if a:
      continue
    else:
      chrom_list.append(chromosome)
  return ' '.join(chrom_list)

def create_chrom_sizes_dict(chrom_sizes_fn):
    ''' it assumes the chrom sizes fn has chr<k>\tsize format 
    and returns a dictionary with key chr<k> and value int(size) '''
    chrom_sizes_dict = {}
    f_in = open(chrom_sizes_fn, 'r')
    for line in f_in:
        cols = line.strip().split('\t')
        chrom = cols[0]
        size = int(cols[1])
        chrom_sizes_dict[chrom] = size
    f_in.close()
    return chrom_sizes_dict

def get_number_of_reads_bam_file(bam_fn):
    ''' count the number of reads in a bam file (fullname) using pysam 
        Note: the bam file needs to be sorted and indexed for this command to 
        work.
    '''
    try:
        n_reads = reduce(lambda x, y: x + y, [ eval('+'.join(l.rstrip('\n').split('\t')[2:]) ) for l in pysam.idxstats(bam_fn) ])
    except:
        n_reads = -1
    return float(n_reads)

def fastx_trimmer(options, input_fn, output_fn):
    ''' writes command for running fastx_trimmer 
    options should be a string -Q 33 -f 16'''
    en = 'fastx_trimmer'
    cmd = '%s %s -i %s -o %s;' %(en, options, input_fn, output_fn)
    return cmd

def casava_quality_filter(input_fn, output_fn):
    ''' writes command for running fastx_trimmer 
    options should be a string -Q 33'''
    en = 'perl /seq/epiprod/de/Cerebellum/code/casava_filter.pl'
    cmd = '%s -i %s -o %s;' %(en, input_fn, output_fn)
    return cmd

def fastx_adaptor_filter(adaptor_seq, options, input_fn, output_fn):
    ''' writes command for running fastx_clipper
    options should be a string -Q 33'''
    en = 'fastx_clipper'
    cmd = '%s -a %s %s -i %s -o %s;' %(en, adaptor_seq, options, input_fn, output_fn)
    return cmd

def fastq_quality_filter(options, input_fn, output_fn):
    ''' writes command for running fastx_trimmer 
    options should be a string -Q 33'''
    en = 'fastq_quality_filter'
    cmd = '%s %s -i %s -o %s;' %(en, options, input_fn, output_fn)
    return cmd

def fastx_artifacts_filter(options, input_fn, output_fn):
    ''' writes command for running fastx_trimmer 
    options should be a string -Q 33'''
    en = 'fastx_artifacts_filter'
    cmd = '%s %s -i %s -o %s;' %(en, options, input_fn, output_fn)
    return cmd

def fastqc(input_fns, output_dir):
    ''' writes command for running fastx_trimmer 
    options should be a string -Q 33'''
    en = 'perl /broad/software/free/Linux/redhat_5_x86_64/pkgs/fastqc_0.10.1/FastQC/fastqc'
    cmd = 'cd %s;%s %s' %(output_dir, en, input_fns)
    return cmd

def trim_galore_filter(adapter_seq, options, input_fn, output_dir):
    ''' writes command for running fastx_trimmer 
    options should be a string -Q 33'''
    en = '/home/unix/dfernand/bin/trim_galore/trim_galore'
    cmd = 'cd %s; %s -a %s %s %s;' %(output_dir, en, adapter_seq, options, input_fn)
    return cmd

def bowtie_1_run(options, index_fn, input_fn, output_fn):
    ''' writes command for running bowtie mapper 
        NOTE: output_fn always needs to be a basefullname with .sorted, it will add the .bam
    '''
    en = "/home/unix/dfernand/bin/bowtie-1.0.0/bowtie"
    cmd = "%s %s %s %s | samtools view -bS - | samtools sort -n - %s" %(en, options, index_fn, input_fn, output_fn)
    return cmd

def list_of_random_integers_no_repetition(x, N):
    ''' returns a list of x random integers from N total integers with no repetition '''
    answer = set()
    sampleSize = x
    answerSize = 0
    while answerSize < sampleSize:
        r = random.randint(0, N)
        if r not in answer:
            answerSize += 1
            answer.add(r)
    return answer 

def bufcount(fullname):
    ''' count number of lines in a file, slower than wccount 
    No good for v large files, too much mem requirement
    https://gist.github.com/zed/0ac760859e614cd03652
    '''
    f = open(fullname)                  
    lines = 0
    buf_size = 1024 * 1024
    read_f = f.read # loop optimization
    buf = read_f(buf_size)
    while buf:
        lines += buf.count('\n')
        buf = read_f(buf_size)
    f.close()
    return lines

def wccount(fullname):
    ''' count number of lines in a file using wc -l 
    https://gist.github.com/zed/0ac760859e614cd03652
    '''
    out = subprocess.Popen(['wc', '-l', fullname],
                         stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT
                         ).communicate()[0]
    return int(out.partition(b' ')[0])

def genome_Nbases(genome_length_fn, assembly):
    ''' file should be chr<i>\tnumber of bases '''
    Nbases = 0
    if assembly == 'hg19':
        chrom_set = set(['chr%s' %(i) for i in range(1,23)]+['chrX','chrY'])
    elif assembly == 'mm9' or assembly == 'mm10':
        chrom_set = set(['chr%s' %(i) for i in range(1,20)]+['chrX','chrY'])
    with open(genome_length_fn, 'r') as f:
        for line in f:
            cols = line.strip().split('\t')
            if cols[0] in chrom_set:
                Nbases += int(cols[1])
    return Nbases

def genome_sizes_dict(genome_length_fn, assembly):
    ''' file should be chr<i>\tnumber of bases '''
    if assembly == 'hg19':
        tmp1 = {'chr'+str(i):0 for i in range(1,23)}
        tmp2 = {'chrX':0, 'chrY':0}
        chrom_sizes_dict = dict(tmp1.items() + tmp2.items())
    elif assembly == 'mm9' or assembly == 'mm10':
        tmp1 = {'chr'+str(i):0 for i in range(1,20)}
        tmp2 = {'chrX':0, 'chrY':0}
        chrom_sizes_dict = dict(tmp1.items() + tmp2.items())
    with open(genome_length_fn, 'r') as f:
        for line in f:
            cols = line.strip().split('\t')
            if cols[0] in chrom_sizes_dict:
                chrom_sizes_dict[cols[0]] = int(cols[1])
    return chrom_sizes_dict

def weighted_sample(items, n):
    ''' weighted sampling with replacement 
    http://stackoverflow.com/questions/2140787/select-random-k-elements-from-a-list-whose-elements-have-weights/2149533#2149533
    Note: it returns a generator object http://stackoverflow.com/questions/102535/what-can-you-use-python-generator-functions-for
    To cath the generator object as a list do a = list(weighted_sample(items, n))
    '''
    total = float(sum(w for w, v in items))
    i = 0
    w, v = items[0]
    while n:
        x = total * (1 - random.random() ** (1.0 / n))
        total -= x
        while x > w:
            x -= w
            i += 1
            w, v = items[i]
        w -= x
        yield v
        n -= 1

#################################################################################
############# ALL THIS IS FOR WEIGTHED SAMPLING WITHOUT REPLACEMENT #############
# http://stackoverflow.com/questions/2140787/select-random-k-elements-from-a-list-whose-elements-have-weights/2149533#2149533
# R function sample(x, size, replace = FALSE, prob = NULL)
#################################################################################
class Node:
    # Each node in the heap has a weight, value, and total weight.
    # The total weight, self.tw, is self.w plus the weight of any children.
    __slots__ = ['w', 'v', 'tw']
    def __init__(self, w, v, tw):
        self.w, self.v, self.tw = w, v, tw

def rws_heap(items):
    # h is the heap. It's like a binary tree that lives in an array.
    # It has a Node for each pair in `items`. h[1] is the root. Each
    # other Node h[i] has a parent at h[i>>1]. Each node has up to 2
    # children, h[i<<1] and h[(i<<1)+1].  To get this nice simple
    # arithmetic, we have to leave h[0] vacant.
    h = [None]                          # leave h[0] vacant
    for w, v in items:
        h.append(Node(w, v, w))
    for i in range(len(h) - 1, 1, -1):  # total up the tws
        h[i>>1].tw += h[i].tw           # add h[i]'s total to its parent
    return h

def rws_heap_pop(h):
    gas = h[1].tw * random.random()     # start with a random amount of gas

    i = 1                     # start driving at the root
    while gas > h[i].w:       # while we have enough gas to get past node i:
        gas -= h[i].w         #   drive past node i
        i <<= 1               #   move to first child
        if gas > h[i].tw:     #   if we have enough gas:
            gas -= h[i].tw    #     drive past first child and descendants
            i += 1            #     move to second child
    w = h[i].w                # out of gas! h[i] is the selected node.
    v = h[i].v

    h[i].w = 0                # make sure this node isn't chosen again
    while i:                  # fix up total weights
        h[i].tw -= w
        i >>= 1
    return v

def random_weighted_sample_no_replacement(items, n):
    heap = rws_heap(items)              # just make a heap...
    for i in range(n):
        yield rws_heap_pop(heap)        # and pop n items off it.


