#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2021 Ryan Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

"""
Create neutral sites mask for athena mutation rate model
"""


import pandas as pd
import pybedtools as pbt
import argparse
from os.path import splitext
from athena.utils.misc import bgzip


def load_constraint(tsv_in):
    """
    Load gnomAD constraint .tsv and return list of non-dispensable transcript IDs
    """

    # Define gnomAD constraint cutoffs
    cutoffs = {
        'oe_mis' : 1,
        'oe_lof' : 1,
        'pLI' : 0.1
    }

    # Load constraint data
    if tsv_in.endswith('gz'):
        df = pd.read_csv(tsv_in, sep='\t', compression='gzip')
    else:
        df = pd.read_csv(tsv_in, sep='\t')

    # Filter constraint data
    keepers = (df.oe_mis < cutoffs['oe_mis']) | \
              (df.oe_lof < cutoffs['oe_mis']) | \
              (df.pLI > cutoffs['pLI'])

    return df.transcript[keepers].tolist()


def get_exons(gtf_in, tx_ids):
    """
    Extract all exons from a GTF that match any transcript in tx_ids
    """

    exons_str = ''

    for feature in pbt.BedTool(gtf_in):
        # Only keep exons
        if feature.fields[2] != 'exon':
            continue

        # Only keep transcripts in tx_ids
        if feature.attrs.get('transcript_id', '').split('.')[0] not in tx_ids:
            continue

        ex_str = '\t'.join(str(x) for x in [feature.chrom, feature.start, feature.end])

        exons_str += ex_str + '\n'

    return pbt.BedTool(exons_str, from_string=True).sort().merge().saveas()


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--gtf', required=True, help='Input GTF')
    parser.add_argument('--constraint', required=True, help='Input gnomAD ' + 
                        'constraint .tsv')
    parser.add_argument('--gerp', required=True, help='Input GERP BED file')
    parser.add_argument('--outfile', required=True, help='Output BED')
    parser.add_argument('-z', '--bgzip', action='store_true', help='Compress ' +
                        'output with bgzip')
    args = parser.parse_args()

    # Load and filter gnomAD constraint data
    tx_ids = load_constraint(args.constraint)

    # Parse gtf and extract exons from transcripts in tx_ids
    exons_bt = get_exons(args.gtf, tx_ids)

    # Load GERP elements
    gerp_bt = pbt.BedTool(args.gerp).sort().merge().saveas()

    # Merge exons & GERP elements
    mask_bt = exons_bt.cat(gerp_bt).sort().merge().saveas()
    
    # Write to file
    outfile = args.outfile
    if outfile.endswith('gz'):
        outfile = splitext(outfile)[0]
    mask_bt.saveas(outfile)

    # Bgzip output, if optioned
    if args.bgzip:
        bgzip(outfile)


if __name__ == '__main__':
    main()

