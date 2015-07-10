# Example usage:
# python gene_analysis.py --gene-file data/flu_gene_ids.tsv Flu --gene-file data/hek_gene_ids.tsv HEK --tsv data/flu-hek_gene_ids.tsv > output/flu-hek_gene_analysis.txt

import click
from collections import defaultdict
from itertools import combinations_with_replacement
from localization import get_lcc_size
from separation import (calc_single_set_distance, calc_set_pair_distances,
    read_gene_list, read_network, remove_self_links)
from StringIO import StringIO

cache = {}
for analysis_file in ('data/flu_gene_analysis.tsv', 'data/flu-hek_gene_analysis.tsv', 'data/flu-jurkat_gene_analysis.tsv', 'data/hek-jurkat_gene_analysis.tsv'):
    s_values = {}
    with open(analysis_file) as f:
        for line in f.read().splitlines():
            if line.startswith('#'):
                continue
            rank, protein_a, protein_b, d_AB, s_AB = line.split('\t')
            cache[(protein_a, protein_b)] = float(d_AB), float(s_AB)

@click.command()
@click.option('--interactome', default='data/interactome.tsv',
              help='Interactome file to use')
@click.option('--tsv', default='',
              help='Outputs a tsv of the results to the given file')
@click.option('--gene-file', default=('data/flu_gene_ids.tsv', 'Flu'),
              multiple=True, type=(str, str))
def main(interactome, tsv, gene_file):
    # Load network
    G = read_network(interactome)
    all_genes_in_network = set(G.nodes())
    remove_self_links(G)

    if tsv:
        with open(tsv, 'wb') as f:
            f.write('#\trank\tprotein_a\tprotein_b\td_AB\ts_AB\n')
    print ''

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

    # Returns d_AB and s_AB of a protein pair
    def analyze_proteins(protein_a, protein_b):
        if (protein_a, protein_b) in cache:
            return cache[(protein_a, protein_b)]

        genes_A = set(genes[protein_a]) & all_genes_in_network
        genes_B = set(genes[protein_b]) & all_genes_in_network

        # Perform calculations
        d_A = calc_single_set_distance(G, genes_A)
        d_B = calc_single_set_distance(G, genes_B)
        d_AB = calc_set_pair_distances(G, genes_A, genes_B)
        s_AB = d_AB - (d_A + d_B)/2.

        return d_AB, s_AB

    # Print information about each protein set as debugging information
    print 'PROTEINS'
    for i, protein in enumerate(sorted(genes)):
        print str(i + 1) + '. ' + protein + '; S = ' + \
            str(get_lcc_size(G, genes[protein]))
    print ''

    # Analyze each protein combination and sort from lowest s_AB to highest
    count = 1
    analyses = []
    for protein_a, protein_b in map(tuple, combinations_with_replacement(sorted(genes), 2)):
        d_AB, s_AB = analyze_proteins(protein_a, protein_b)
        analyses.append((protein_a, protein_b, d_AB, s_AB))
        count += 1
    analyses.sort(key=lambda x: x[3])

    # Print analyses in order
    for i, analysis in enumerate(analyses):
        if tsv:
            with open(tsv, 'ab') as f:
                f.write('\t'.join((str(i + 1), analysis[0], analysis[1],
                        str(analysis[2]), str(analysis[3]))) + '\n')
        print str(i + 1) + '. Proteins: ' + analysis[0] + ' and ' \
            + analysis[1]
        print 'd_AB = ' + str(analysis[2])
        print 's_AB = ' + str(analysis[3])
        print ''

if __name__ == '__main__':
    main()
