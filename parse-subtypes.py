from Bio import SeqIO
import argparse
import re

pat = re.compile('(genotype|subtype)[:]*[ ]*([1-7][a-z])')

parser = argparse.ArgumentParser(
    description='Extract HCV genotype/subtype information.'
)
parser.add_argument('gb', help='Path to GenBank formatted file.')
parser.add_argument('csv', help='Path to write CSV output.')

args = parser.parse_args()

outfile = open(args.csv, 'w')
outfile.write('accession,gb_genotype\n')

gb = SeqIO.parse(args.gb, 'genbank')
for i, record in enumerate(gb):
    # attempt to extract genotype from source annotation
    source = None
    for feature in record.features:
        if feature.type != 'source':
            continue
        source = feature

    note = None
    if source:
        note = source.qualifiers.get('note', [None])[0]

    genotype = ''
    if note:
        match = pat.findall(note)
        if len(match) > 0:
            genotype = match[0][1]

    # attempt to extract genotype from record description
    if len(genotype) == 0:
        match = pat.findall(record.description)
        if len(match) > 0:
            genotype = match[0][1]

    outfile.write('%s,%s\n' % (record.id, genotype))

outfile.close()
