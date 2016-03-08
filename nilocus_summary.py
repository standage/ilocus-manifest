#!/usr/bin/env python

from __future__ import division
from __future__ import print_function
import argparse
import pandas
import re
import sys

def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--gff3', type=argparse.FileType('r'),
                        help='validate total genome size with iLocus GFF3')
    parser.add_argument('iloc', type=argparse.FileType('r'),
                        help='iLocus table')
    return parser


def calc_genome_size(gff3):
    genomesize = 0
    for line in gff3:
        if not line.startswith('##sequence-region'):
            continue

        match = re.search('##sequence-region\s+(\S+)\s+(\d+)\s+(\d+)', line)
        errmsg = 'unable to parse sequence-region pragma: {}'.format(line)
        assert match, errmsg

        seqid = match.group(1)
        start = int(match.group(2))
        end = int(match.group(3))
        errmsg = ('expected full genomic sequences, but sequence {} '
                  'starts with nucleotide {} != 1'.format(seqid, start))
        assert start == 1, errmsg
        genomesize += end

    return genomesize


def main(args):
    gff3_genomesize = 0
    if args.gff3:
        gff3_genomesize = calc_genome_size(args.gff3)
    iloci = pandas.read_table(args.iloc)
    table_genomesize = iloci['EffectiveLength'].sum()
    if gff3_genomesize:
        errmsg = ('genome size mismatch: {} != {}'.format(gff3_genomesize,
                                                          table_genomesize))
        assert gff3_genomesize == table_genomesize, errmsg

    niloci = iloci.loc[(iloci.LocusClass == 'niLocus')]
    if len(niloci) == 0:
        print('No niLoci!')
        sys.exit(0)
    print('niLoci', 'SpaceOccupied', 'Length', '%GC Content', sep='\t')
    spaceocc = niloci['EffectiveLength'].sum()
    spaceoccperc = spaceocc / table_genomesize * 100
    occ = '{:.1f} Mb ({:.1f}%)'.format(spaceocc / 1000000, spaceoccperc)

    lenqnt = list(niloci['Length'].quantile([.25, .5, .75]))
    lenqnt = [int(x) for x in lenqnt]
    gcqnt = list(niloci['GCContent'].quantile([.25, .5, .75]))
    gcqnt = [float('{:.3f}'.format(x)) for x in gcqnt]
    print(len(niloci), occ, lenqnt, gcqnt, sep='\t')


if __name__ == '__main__':
    main(get_parser().parse_args())


