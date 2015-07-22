import math
import re

def num_genes(protein):
    if protein.endswith('[Flu]'):
        with open('data/flu_gene_ids.tsv') as f:
            return len(filter(lambda x: not x.startswith('#')
                                  and x.split('\t')[0] + ' [Flu]' == protein,
                              f.read().splitlines()))
    elif protein.endswith(' [Jurkat]'):
        with open('data/jurkat_gene_ids.tsv') as f:
            return len(filter(lambda x: not x.startswith('#')
                                  and x.split('\t')[0] + ' [Jurkat]' == protein,
                              f.read().splitlines()))

with open('data/flu-jurkat_gene_analysis.tsv') as f:
    s_values = {}
    for line in f.read().splitlines():
        if line.startswith('#'):
            continue
        rank, protein_a, protein_b, d_AB, s_AB = line.split('\t')
        s_values[tuple(sorted((protein_a, protein_b)))] = float(s_AB)

    normalized_values = {}
    for key in s_values:
        protein_a, protein_b = key
        normalized_values[key] = s_values[key] / math.sqrt(num_genes(protein_a) * num_genes(protein_b))

    print '\t'.join(('Rank', 'Flu Protein', 'Jurkat Protein', 'Size-normalized s_AB', 's_AB'))
    i = 1
    for protein_a, protein_b in sorted(normalized_values, key=normalized_values.get):
        disease_a = re.match(r'.+ \[(.*)\]', protein_a).group(1)
        disease_b = re.match(r'.+ \[(.*)\]', protein_b).group(1)
        if disease_a != disease_b:
            if protein_a.endswith('[Flu]'):
                protein_flu, protein_jurkat = protein_a[:-6], protein_b[:-9]
            else:
                protein_flu, protein_jurkat = protein_b[:-6], protein_a[:-9]
            print '\t'.join((str(i), protein_flu, protein_jurkat, str(normalized_values[(protein_a, protein_b)]), str(s_values[(protein_a, protein_b)])))
            i += 1
