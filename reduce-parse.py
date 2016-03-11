"""
Exclude clipped gene sequences from parse-genbank.csv output.
"""

from csv import DictReader, DictWriter

reader = DictReader(open('/Users/art/wip/HCVdb/parse-genbank.csv', 'rU'))
outfile = DictWriter(open('/Users/art/wip/HCVdb/HCV_AnyRegion.annotations.csv', 'w'),
                     fieldnames=['accession', 'description', 'coldate', 'country', 'note'],
                     extrasaction='ignore'
                     )
outfile.writeheader()
for row in reader:
    outfile.writerow(row)

