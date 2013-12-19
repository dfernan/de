#!/usr/bin/env python

# Import modules
import optparse

def split_vcf(vcf_fn, vcf_c57_out_fn, vcf_cast_out_fn):
    '''
    This code takes a diploid VCF file and splits it into two haploid 
    VCFs
    ONLY optimized for use within 
    '''
    vcf_in = open(vcf_fn,'r')
    vcf_c57_out = open(vcf_c57_out_fn,'w')
    vcf_cast_out = open(vcf_cast_out_fn,'w')
    for line in vcf_in:
        cols = line.split('\t')
        geno_index = cols[9]
        geno_c57_index = int(geno_index[0])
        geno_cast_index = int(geno_index[2])
        nt = [cols[3]]+cols[4].split(',')
        if geno_c57_index != 0:
            cols_c57 = cols
            cols_c57[4] = nt[geno_c57_index]
            new_geno = list(cols_c57[9])
            new_geno[0] = '1'
            new_geno[2] = '1'
            cols_c57[9] = ''.join(new_geno)
            vcf_c57_out.write('\t'.join(cols_c57))
        if geno_cast_index != 0:
            cols_cast = cols
            cols_cast[4] = nt[geno_cast_index]
            new_geno = list(cols_cast[9])
            new_geno[0] = '1'
            new_geno[2] = '1'
            cols_cast[9] = ''.join(new_geno)
            vcf_cast_out.write('\t'.join(cols_cast))
    vcf_in.close()
    vcf_c57_out.close()
    vcf_cast_out.close()
    return 0

def main():
    parser = optparse.OptionParser()
    parser.add_option('-i', '--vcf_fn', action = "store")
    parser.add_option('-c', '--vcf_c57_out_fn', action = "store")
    parser.add_option('-a', '--vcf_cast_out_fn', action = "store")
    myargs, noidea = parser.parse_args()
    split_vcf(myargs.vcf_fn, myargs.vcf_c57_out_fn, myargs.vcf_cast_out_fn)
    return 0

if __name__ == "__main__":
    main()
