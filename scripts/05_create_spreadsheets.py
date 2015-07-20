from collections import namedtuple
import pandas as pd
import re

def parse_protein(protein):
    m = re.search('(?P<protein_name>.+) \[(?P<disease>.*)\]', protein)
    return {'protein': m.group('protein_name'), 'disease': m.group('disease')}

# Non-normalized data
with open('data/flu-jurkat_crossed_gene_analysis.tsv') as f:
    with open('output/flu-jurkat_crossed.tsv', 'wb') as out:
        out.write('\t'.join(('Rank', 'Flu Protein', 'Jurkat Protein', 's_AB')) + '\n')
        for line in f.read().splitlines():
            if line.startswith('#'):
                continue
            rank, protein_a, protein_b, d_AB, s_AB = line.split('\t')
            protein_a = parse_protein(protein_a)
            protein_b = parse_protein(protein_b)
            if protein_a['disease'] == 'Jurkat':
                protein_a, protein_b = protein_b, protein_a
            assert protein_a['disease'] == 'Flu' and protein_b['disease'] == 'Jurkat'
            out.write('\t'.join((rank, protein_a['protein'], protein_b['protein'], s_AB)))
            out.write('\n')

# Normalized data
data = pd.DataFrame.from_csv('output/flu-jurkat_normalized_crossed_matrix.tsv', sep='\t').to_dict()
s_values = []
for flu_protein, s_data in data.iteritems():
    for jurkat_protein, s_AB in s_data.iteritems():
        flu = parse_protein(flu_protein)['protein']
        jurkat = parse_protein(jurkat_protein)['protein']
        s_values.append((flu, jurkat, float(s_AB)))
with open('output/flu-jurkat_normalized_crossed.tsv', 'wb') as out:
    out.write('\t'.join(('Rank', 'Flu Protein', 'Jurkat Protein', 'Normalized s_AB')) + '\n')
    for i, data_row in enumerate(sorted(s_values, key=lambda x: x[2], reverse=True)):
        flu_protein, jurkat_protein, s_AB = data_row
        out.write('\t'.join((str(i + 1), flu_protein, jurkat_protein, str(s_AB))) + '\n')
