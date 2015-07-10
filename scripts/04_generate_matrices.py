from collections import defaultdict
import numpy as np
import pandas as pd
import re

def generate_matrix(input, output, crossed, normalize):
    genes = defaultdict(set)
    s_values = {}
    with open(input) as f:
        for line in f.read().splitlines():
            if line.startswith('#'):
                continue
            rank, protein_a, protein_b, d_AB, s_AB = line.split('\t')
            genes[re.match(r'.+ \[(.*)\]', protein_a).group(1)].add(protein_a)
            genes[re.match(r'.+ \[(.*)\]', protein_b).group(1)].add(protein_b)
            s_values[tuple(sorted((protein_a, protein_b)))] = float(s_AB)
    assert len(genes) == 2
    for protein in genes:
        genes[protein] = list(sorted(genes[protein]))

    if normalize:
        normalized_values = {}
        for protein_a, protein_b in s_values:
            key = tuple(sorted((protein_a, protein_b)))
            normalized_values[key] = s_values[key] / max(s_values[(protein_a, protein_a)], s_values[(protein_b, protein_b)])

    if crossed:
        column_labels = sorted(genes[list(genes.viewkeys())[0]])
        row_labels = sorted(genes[list(genes.viewkeys())[1]])
    else:
        column_labels = row_labels = sorted(genes[list(genes.viewkeys())[0]] + genes[list(genes.viewkeys())[1]])
    data = np.zeros((len(row_labels), len(column_labels)))
    for i, protein_a in enumerate(row_labels):
        for j, protein_b in enumerate(column_labels):
            values = normalized_values if normalize else s_values
            data[i,j] = values[tuple(sorted((protein_a, protein_b)))]

    df = pd.DataFrame(data, index=row_labels, columns=column_labels)
    df.to_csv(output, index=True, header=True, sep='\t')

generate_matrix('data/flu-hek_gene_analysis.tsv', 'output/flu-hek_crossed_matrix.tsv', True, False)
generate_matrix('data/flu-jurkat_gene_analysis.tsv', 'output/flu-jurkat_crossed_matrix.tsv', True, False)
generate_matrix('data/hek-jurkat_gene_analysis.tsv', 'output/hek-jurkat_crossed_matrix.tsv', True, False)
generate_matrix('data/flu-hek_gene_analysis.tsv', 'output/flu-hek_matrix.tsv', False, False)
generate_matrix('data/flu-jurkat_gene_analysis.tsv', 'output/flu-jurkat_matrix.tsv', False, False)
generate_matrix('data/hek-jurkat_gene_analysis.tsv', 'output/hek-jurkat_matrix.tsv', False, False)

generate_matrix('data/flu-hek_gene_analysis.tsv', 'output/flu-hek_normalized_crossed_matrix.tsv', True, True)
generate_matrix('data/flu-jurkat_gene_analysis.tsv', 'output/flu-jurkat_normalized_crossed_matrix.tsv', True, True)
generate_matrix('data/hek-jurkat_gene_analysis.tsv', 'output/hek-jurkat_normalized_crossed_matrix.tsv', True, True)
generate_matrix('data/flu-hek_gene_analysis.tsv', 'output/flu-hek_normalized_matrix.tsv', False, True)
generate_matrix('data/flu-jurkat_gene_analysis.tsv', 'output/flu-jurkat_normalized_matrix.tsv', False, True)
generate_matrix('data/hek-jurkat_gene_analysis.tsv', 'output/hek-jurkat_normalized_matrix.tsv', False, True)
