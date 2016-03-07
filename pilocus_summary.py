#!/usr/bin/env python

from __future__ import division
from __future__ import print_function
import argparse
import pandas
import re

def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--gff3', type=argparse.FileType('r'),
                        help='validate total genome size with iLocus GFF3')
    parser.add_argument('iloc', type=argparse.FileType('r'),
                        help='iLocus table')
    parser.add_argument('premrna', type=argparse.FileType('r'),
                        help='pre-mRNA table')
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

    siloci = iloci.loc[(iloci.LocusClass == 'piLocus')]
    ciloci = iloci.loc[(iloci.LocusClass == 'ciLocus')]
    piloci = pandas.concat([siloci, ciloci])
    premrnas = pandas.read_table(args.premrna)

    print('siLoci', 'ciLoci', 'SpaceOccupied', 'Length', 'ExonCount', '%GC Content', sep='\t')
    spaceocc = piloci['EffectiveLength'].sum()
    spaceoccperc = spaceocc / table_genomesize * 100
    occ = '{:.1f} Mb ({:.1f}%)'.format(spaceocc / 1000000, spaceoccperc)

    lenqnt = list(piloci['Length'].quantile([.25, .5, .75]))
    lenqnt = [int(x) for x in lenqnt]
    exonqnt = list(premrnas['ExonCount'].quantile([.25, .5, .75]))
    exonqnt = [int(x) for x in exonqnt]
    gcqnt = list(piloci['GCContent'].quantile([.25, .5, .75]))
    gcqnt = [float('{:.3f}'.format(x)) for x in gcqnt]
    print(len(siloci), len(ciloci), occ, lenqnt, exonqnt, gcqnt, sep='\t')


if __name__ == '__main__':
    main(get_parser().parse_args())


