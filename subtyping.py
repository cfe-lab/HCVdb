from Bio import SeqIO
import argparse
import subprocess

import re
cigar_re = re.compile('[0-9]+[MIDNSHPX=]')  # CIGAR token

def match_length(cigar):
    """Returns total number of bases matched to a reference."""
    sizes = re.findall(r'(\d+)M', cigar)
    return sum(map(int, sizes))


parser = argparse.ArgumentParser(
    description='Use bowtie2 to map sequences to HCV reference genotypes '
                'and subtypes.'
)
parser.add_argument('infile', help='Path to FASTA input.')
parser.add_argument('index', help='Path to Bowtie2 index files.')
parser.add_argument('outfile', help='Path to write CSV output.')

args = parser.parse_args()

p = subprocess.Popen([
    'bowtie2',
    '-f', args.infile,
    '-x', args.index,
    '--local',
    '--no-hd',  # no header lines
    '--quiet',
    '-p6'  # number of threads
], stdout=subprocess.PIPE)

outfile = open(args.outfile, 'w')
outfile.write('qname,flag,rname,pos,mapq,match.len,seq.len\n')

for line in p.stdout:
    qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual = line.strip('\n').split('\t')[:11]
    seqlen = len(seq)
    subtype = ''
    if rname != '*':
        subtype = rname.split('.')[1]

    outfile.write(','.join(map(str, [qname, flag, subtype, pos, mapq, match_length(cigar), seqlen])))
    outfile.write('\n')

p.wait()
if p.returncode():
    raise subprocess.CalledProcessError(p.returncode, 'bowtie2')

outfile.close()
