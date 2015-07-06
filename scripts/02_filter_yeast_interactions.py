# Filters the interactome.tsv to include only yeast two-hybrid interactions,
# and creates a new tsv from them
#
# Usage: python scripts/02_filter_yeast_interactions > data/y2h-interactome.tsv

import sys

print '# This file contains a filtered version of the Human Interactome that includes'
print '# only the yeast two-hybrd interactions ("Binary interactions").'
print '# gene_ID_1	 gene_ID_2	data_source(s)'

with open('data/interactome.tsv') as f:
    for line in f:
        # Ignore comments
        if line.startswith('#'):
            continue
        if 'binary' in line.split('\t')[2].split(';'):
            sys.stdout.write(line)
