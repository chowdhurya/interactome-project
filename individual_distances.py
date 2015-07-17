# Example usage:
# python gene_analysis.py --gene-file data/flu_gene_ids.tsv Flu --gene-file data/hek_gene_ids.tsv HEK --tsv data/flu-hek_gene_ids.tsv > output/flu-hek_gene_analysis.txt

import click
from collections import defaultdict
from itertools import combinations_with_replacement
from localization import get_lcc_size
import numpy as np
import re
from separation import (calc_single_set_distance, calc_set_pair_distances,
    read_gene_list, read_network, remove_self_links, get_pathlengths_for_two_sets)
from StringIO import StringIO
import sys

def _title_from_protein(protein):
    return re.match(r'.+ \[(.*)\]', protein).group(1)

@click.command()
@click.option('--interactome', default='data/interactome.tsv',
              help='Interactome file to use')
@click.option('--gene-file', default=('', 'Flu'), multiple=True,
              type=(str, str))
def main(interactome, gene_file):
    # Load network
    sys.stdout = StringIO()
    G = read_network(interactome)
    all_genes_in_network = set(G.nodes())
    remove_self_links(G)
    sys.stdout = sys.__stdout__

    # Read gene files and separate by protein
    genes = defaultdict(list)
    for data_file, title in gene_file:
        with open(data_file) as f:
            for line in f.read().splitlines():
                if line.startswith('#'):
                    continue
                protein, gene_ids = line.split('\t')
                if title and len(gene_file) > 1:
                    protein = '%s [%s]' % (protein, title)
                genes[protein].extend(gene_ids.split(';'))

    # Returns d_AB, s_AB, and ((gene_1, gene_2, dist)...) of a protein pair
    def analyze_proteins(protein_a, protein_b):
        genes_A = set(genes[protein_a]) & all_genes_in_network
        genes_B = set(genes[protein_b]) & all_genes_in_network

        all_path_lengths = get_pathlengths_for_two_sets(G, genes_A, genes_B)
        all_distances = []

        # Perform calculations
        d_A = calc_single_set_distance(G, genes_A)
        d_B = calc_single_set_distance(G, genes_B)

        for gene_A in genes_A:
            all_distances_A = []
            for gene_B in genes_B:
                if gene_A == gene_B:
                    all_distances_A.append((gene_A, gene_B, 0))
                else:
                    try:
                        all_distances_A.append((gene_A, gene_B, all_path_lengths[min(gene_A, gene_B)][max(gene_A, gene_B)]))
                    except KeyError:
                        pass
            if len(all_distances_A) > 0:
                all_distances.append(min(all_distances_A, key=lambda x: x[2]))

        for gene_B in genes_B:
            all_distances_B = []
            for gene_A in genes_A:
                if gene_A == gene_B:
                    all_distances_B.append((gene_A, gene_B, 0))
                else:
                    try:
                        all_distances_B.append((gene_B, gene_A, all_path_lengths[min(gene_A, gene_B)][max(gene_A, gene_B)]))
                    except KeyError:
                        pass
            if len(all_distances_B) > 0:
                all_distances.append(min(all_distances_B, key=lambda x: x[2]))

        d_AB = np.mean(map(lambda x: x[2], all_distances))
        s_AB = d_AB - (d_A + d_B)/2.

        return d_AB, s_AB, sorted(all_distances, key=lambda x: x[2])

    # Print information about each protein set as debugging information
    print 'PROTEINS'
    for i, protein in enumerate(sorted(genes)):
        print str(i + 1) + '. ' + protein + '; S = ' + \
            str(get_lcc_size(G, genes[protein]))
    print ''

    # Analyze each protein combination and sort from lowest s_AB to highest
    analyses = []
    for protein_a, protein_b in map(tuple, combinations_with_replacement(sorted(genes), 2)):
        d_AB, s_AB, distances = analyze_proteins(protein_a, protein_b)
        analyses.append((protein_a, protein_b, d_AB, s_AB, distances))
    analyses.sort(key=lambda x: x[3])

    # Print analyses in order
    for i, analysis in enumerate(analyses):
        protein_a, protein_b, d_AB, s_AB, distances = analysis
        print str(i + 1) + '. Proteins: ' + protein_a + ' and ' + protein_b
        print 'd_AB = ' + str(d_AB)
        print 's_AB = ' + str(s_AB)
        for gene_1, gene_2, dist in distances:
            print gene_1 + ' -> ' + gene_2 + ' = ' + str(dist)
        print ''

if __name__ == '__main__':
    main()
