#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2021 Ryan Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

"""
Clean all allele frequency annotations from an input VCF
"""


import sys
import argparse
import pysam


# Define INFO entries to clean
freq_infos = 'AN AC AF N_BI_GENOS N_HOMREF N_HET N_HOMALT FREQ_HOMREF ' + \
             'FREQ_HET FREQ_HOMALT CN_NUMBER CN_COUNT CN_FREQ ' + \
             'CN_NONREF_COUNT CN_NONREF_FREQ'
freq_infos = freq_infos.split()


def clean_af_info(record):
    """
    Remove all allele frequency info from a single record
    """

    fields_to_pop = set()
    for field in record.info.keys():
        if field in freq_infos \
        or len([f for f in freq_infos if '_' + f in field]) > 0:
            fields_to_pop.add(field)

    for field in fields_to_pop:
        record.info.pop(field)

    return record


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf', help='Input vcf. Also accepts "stdin" and "-".')
    parser.add_argument('fout', help='Output vcf. Also accepts "stdout" and "-".')
    parser.add_argument('-v', '--verbose', action='store_true', 
                        help='Print progress.')
    args = parser.parse_args()

    # Open connection to input VCF
    if args.vcf in '- stdin'.split():
        vcf = pysam.VariantFile(sys.stdin) 
    else:
        vcf = pysam.VariantFile(args.vcf)
    
    # Open connection to output VCF
    if args.fout in '- stdout'.split():
        fout = pysam.VariantFile(sys.stdout, 'w', header=vcf.header)
    else:
        fout = pysam.VariantFile(args.fout, 'w', header=vcf.header)

    # Clean each record before writing to file
    if args.verbose:
        k, n = 0, 0
    for record in vcf.fetch():
        cleaned_record = clean_af_info(record)
        fout.write(cleaned_record)
        if args.verbose:
            k += 1; n += 1
            if k == 1000 and args.verbose:
                print('Cleaned {:,} records...'.format(n))

    fout.close()


if __name__ == '__main__':
    main()

