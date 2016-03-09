from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(
    description='Convert GenBank records into FASTA format.'
)
parser.add_argument('infile', help='Path to GenBank formatted file.')
parser.add_argument('outfile', help='Path to write FASTA output.')

args = parser.parse_args()

gb = SeqIO.parse(args.infile, 'genbank')

with open(args.outfile, 'w') as of:
    for record in gb:
        accno = record.id
        seq = str(record.seq)
        of.write('>%s\n%s\n' % (accno, seq))

# bowtie2 -f HCV_AnyRegion.fa -p6 -x /Users/art/git/HCVdb/data/lanl -S HCV_AnyRegion.sam --no-head --local

# LANL file contains one blank "Ref." sequence label; this corresponds to QC103, an HCV-1d strain
