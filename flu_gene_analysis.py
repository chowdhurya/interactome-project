# python flu_gene_analysis.py > output/flu_gene_analysis.txt

from collections import defaultdict
from itertools import combinations
from localization import get_lcc_size
from separation import (calc_single_set_distance, calc_set_pair_distances,
    read_gene_list, read_network, remove_self_links)

# Load network
G = read_network('data/interactome.tsv')
all_genes_in_network = set(G.nodes())
remove_self_links(G)
print ''

# Read flu genes file and separate by protein
genes = defaultdict(list)
with open('data/flu_gene_ids.tsv') as f:
    for line in f.readlines()[2:]:
        protein, gene_ids = line.split('\t')
        genes[protein].extend(map(lambda x: x.strip(), gene_ids.split(';')))

# Returns d_AB and s_AB of a protein pair
def analyze_proteins(protein_a, protein_b):
    genes_A = set(genes[protein_a]) & all_genes_in_network
    genes_B = set(genes[protein_b]) & all_genes_in_network

    # Perform calculations
    d_A = calc_single_set_distance(G, genes_A)
    d_B = calc_single_set_distance(G, genes_B)
    d_AB = calc_set_pair_distances(G, genes_A, genes_B)
    s_AB = d_AB - (d_A + d_B)/2.

    return d_AB, s_AB

# Returns S, d_s of a single protein
def analyze_single_protein(protein):
    return

# Print information about each protein set
print 'PROTEINS'
for i, protein in enumerate(sorted(genes)):
    print str(i + 1) + '. ' + protein + '; S = ' + \
        str(get_lcc_size(G, genes[protein]))
print ''

# Analyze each protein combination and sort from lowest s_AB to highest
analyses = []
for protein_a, protein_b in map(tuple, combinations(sorted(genes), 2)):
    d_AB, s_AB = analyze_proteins(protein_a, protein_b)
    analyses.append((protein_a, protein_b, d_AB, s_AB))
analyses.sort(key=lambda x: x[3])

# Print analyses in order
for i, analysis in enumerate(analyses):
    print str(i + 1) + '. Proteins: ' + analysis[0] + ' and ' + analysis[1]
    print 'd_AB = ' + str(analysis[2])
    print 's_AB = ' + str(analysis[3])
    print ''
