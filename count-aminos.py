import HyGo
import argparse
from csv import DictReader
from translation import translate

parser = argparse.ArgumentParser(
    description='Count amino acid frequencies from CSV output of parse-genbank.py. '
                'Record insertions relative to H77 in a separate column.'
)
parser.add_argument('infile', help='Path to input CSV generated by parse-genbank.py.')
parser.add_argument('outfile', help='Path to write CSV output.')
args = parser.parse_args()

hygo = HyGo.HyGo()

# change settings for protein alignment
hygo.set_alphabet(HyGo.amino_alphabet)
hygo.set_matrix(HyGo.blosum62)
hygo.set_gap_open(7, True)
hygo.set_gap_open(7, False)
hygo.set_gap_extend(1, True)
hygo.set_gap_extend(1, False)

h77 = {}
with open('data/h77-genes.csv', 'rU') as f:
    for line in f:
        gene, seq = line.strip('\n').split(',')
        h77.update({gene: translate(seq, 0)})

reader = DictReader(open(args.infile, 'rU'))

for row in reader:
    accno = row['accession']
    for gene, refseq in h77.iteritems():
        seq = row[gene]
        if len(seq) == 0:
            # skip genes with no coverage in record
            continue

        best_frame = 0
        best_results = {'score': -1000}
        for frame in range(3):
            prot = translate(seq, frame)
            results = hygo.align(refseq, prot)
            if results['score'] > best_results['score']:
                best_results = results
                best_frame = frame

        # TODO: process insertions

        # TODO: detect frame-shifts


    break
