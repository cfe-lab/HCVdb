from Bio import SeqIO
import HyGo
import argparse
import sys
import os
from csv import DictWriter

hygo = HyGo.HyGo()


parser = argparse.ArgumentParser(
    description='Read records from GenBank file.  Align sequences against HXB2 gene reference '
                'sequences and write to separate files.'
)
parser.add_argument('gb', help='Path to GenBank formatted file.')
parser.add_argument('csv', help='Path to write CSV output.')
parser.add_argument('-cutoff', help='Normalized alignment score cutoff.  Defaults to 1.0.')

args = parser.parse_args()


h77 = {}
with open('data/h77-genes.csv', 'rU') as f:
    for line in f:
        gene, seq = line.strip('\n').split(',')
        h77.update({gene: seq})

gap_open = 15
gap_extend = 5
use_terminal_gap_penalty = 1

# set up output
genes = h77.keys()
genes.sort()
fieldnames = ['accession', 'description', 'coldate', 'country', 'note'] + genes

writer = DictWriter(open(args.csv, 'w'), fieldnames=fieldnames)
writer.writeheader()

gb = SeqIO.parse(args.gb, 'genbank')
for i, record in enumerate(gb):
    if i % 10 == 0:
        print i

    desc = record.description
    if 'Patent' in desc or desc.startswith('Modified') or 'IRNA' in desc.upper() or 'UTR' in desc.split():
        # ignore records associated with a patent
        continue

    seq = str(record.seq)

    source = None
    for feature in record.features:
        if feature.type != 'source':
            continue
        source = feature

    segments = {}
    for gene, refseq in h77.iteritems():
        results = hygo.align(refseq, seq)
        aref = results['ref']
        aquery = results['query']
        ascore = results['score']

        left, right = HyGo.get_boundaries(aref)

        # diagnostics
        overlap = sum(map(lambda i: aquery[i]!='-' and aref[i]!='-', range(left, right)))
        if overlap < 60:
            segments.update({gene: ''})
            continue

        relscore = ascore / float(overlap)
        if relscore < 0.8:
            segments.update({gene: ''})
            continue

        #print ','.join(map(str, [record.id, gene, overlap, relscore]))
        segments.update({gene: aquery[left:right].strip('-')})

    row = {
        'accession': record.id,
        'description': desc,
        'coldate': source.qualifiers.get('collection_date', [''])[0] if source else '',
        'country': source.qualifiers.get('country', [''])[0] if source else '',
        'note': source.qualifiers.get('note', [''])[0] if source else ''
    }
    row.update(segments)
    writer.writerow(row)

