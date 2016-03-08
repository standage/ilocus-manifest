#!/usr/bin/env python

from __future__ import division
from __future__ import print_function
import argparse
import pandas
import re

def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--gff3', type=argparse.FileType('r'), metavar='FILE',
                        help='validate total genome size with iLocus GFF3')
    parser.add_argument('--premrnas', type=argparse.FileType('r'),
                        metavar='FILE',
                        help='pre-mRNA data table (required for piLocus stats')
    parser.add_argument('--type', default='piLocus',
                        choices=['piLocus', 'niLocus', 'iiLocus', 'miLocus'],
                        help='the type of iLoci to summarize')
    parser.add_argument('iloci', type=argparse.FileType('r'),
                        help='iLocus/miLocus table')
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


def summarize_piloci(ilocus_infile, premrna_infile, genomesize=0):
    iloci = pandas.read_table(ilocus_infile)
    if genomesize != 0:
        testsize = iloci['EffectiveLength'].sum()
        errmsg = ('genome size mismatch: {} != {}'.format(genomesize,
                                                          testsize))
        assert genomesize == testsize, errmsg

    print('siLoci', 'ciLoci', 'SpaceOccupied', 'Length', 'ExonCount',
          '%GC Content', sep='\t')

    siloci = iloci.loc[(iloci.LocusClass == 'piLocus')]
    ciloci = iloci.loc[(iloci.LocusClass == 'ciLocus')]
    piloci = pandas.concat([siloci, ciloci])
    spaceocc = piloci['EffectiveLength'].sum()
    
    premrnas = pandas.read_table(premrna_infile)
    spaceoccperc = spaceocc / genomesize * 100
    occ = '{:.1f} Mb ({:.1f}%)'.format(spaceocc / 1000000, spaceoccperc)

    lenqnt = list(piloci['Length'].quantile([.25, .5, .75]))
    lenqnt = [int(x) for x in lenqnt]
    exonqnt = list(premrnas['ExonCount'].quantile([.25, .5, .75]))
    exonqnt = [int(x) for x in exonqnt]
    gcqnt = list(piloci['GCContent'].quantile([.25, .5, .75]))
    gcqnt = [float('{:.3f}'.format(x)) for x in gcqnt]
    print(len(siloci), len(ciloci), occ, lenqnt, exonqnt, gcqnt, sep='\t')


def summarize_miloci(ilocus_infile, genomesize=0):
    iloci = pandas.read_table(ilocus_infile)
    if genomesize != 0:
        testsize = iloci['EffectiveLength'].sum()
        errmsg = ('genome size mismatch: {} != {}'.format(genomesize,
                                                          testsize))
        assert genomesize == testsize, errmsg

    miloci = iloci.loc[(iloci.LocusClass == 'miLocus')]
    giloci = iloci.loc[(iloci.LocusClass.isin(['piLocus', 'ciLocus', 'niLocus']))]
    print('miLoci', 'SpaceOccupied', 'Length', '%GC Content', 'giLociSingletons', sep='\t')

    spaceocc = miloci['EffectiveLength'].sum()
    spaceoccperc = spaceocc / genomesize * 100
    occ = '{:.1f} Mb ({:.1f}%)'.format(spaceocc / 1000000, spaceoccperc)

    lenqnt = list(miloci['Length'].quantile([.25, .5, .75]))
    lenqnt = [int(x) for x in lenqnt]
    gcqnt = list(miloci['GCContent'].quantile([.25, .5, .75]))
    gcqnt = [float('{:.3f}'.format(x)) for x in gcqnt]
    print(len(miloci), occ, lenqnt, gcqnt, len(giloci), sep='\t')


def summarize_iiloci_niloci(ilocus_infile, genomesize=0, prefix='i'):
    iloci = pandas.read_table(ilocus_infile)
    if genomesize != 0:
        testsize = iloci['EffectiveLength'].sum()
        errmsg = ('genome size mismatch: {} != {}'.format(genomesize,
                                                          testsize))
        assert genomesize == testsize, errmsg

    ilocus_type = prefix + 'iLocus'
    iloci = iloci.loc[(iloci.LocusClass == ilocus_type)]
    print(prefix + 'iLoci', 'SpaceOccupied', 'Length', '%GC Content', sep='\t')
    if len(iloci) == 0:
        print('0', '0.0 Mb (0%)', '-', '-', sep='\t')
        return

    spaceocc = iloci['EffectiveLength'].sum()
    spaceoccperc = spaceocc / genomesize * 100
    occ = '{:.1f} Mb ({:.1f}%)'.format(spaceocc / 1000000, spaceoccperc)

    lenqnt = list(iloci['Length'].quantile([.25, .5, .75]))
    lenqnt = [int(x) for x in lenqnt]
    gcqnt = list(iloci['GCContent'].quantile([.25, .5, .75]))
    gcqnt = [float('{:.3f}'.format(x)) for x in gcqnt]
    print(len(iloci), occ, lenqnt, gcqnt, sep='\t')


def main(args):
    genomesize = 0
    if args.gff3:
        genomesize = calc_genome_size(args.gff3)
    if args.type == 'piLocus':
        summarize_piloci(args.iloci, args.premrnas, genomesize)
    elif args.type == 'niLocus':
        summarize_iiloci_niloci(args.iloci, genomesize, prefix='n')
    elif args.type == 'iiLocus':
        summarize_iiloci_niloci(args.iloci, genomesize, prefix='i')
    elif args.type == 'miLocus':
        summarize_miloci(args.iloci, genomesize)

if __name__ == '__main__':
    main(get_parser().parse_args())


