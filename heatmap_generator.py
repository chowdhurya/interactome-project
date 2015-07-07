import matplotlib
matplotlib.use('Agg')
import click
import matplotlib.pyplot as plt
import numpy as np

@click.command()
@click.option('--input', default='data/flu_gene_analysis.tsv',
              help='Input gene analysis file')
@click.option('--output', default='output/flu_gene_heatmap.png',
              help='Output image destination')
@click.option('--title', default='Flu Gene Analysis', help='Heatmap title')
def main(input, output, title):
    proteins = set()
    s_values = {}
    with open(input) as f:
        for line in f:
            if line.startswith('#'):
                continue
            values = line.split('\t')
            protein_a, protein_b, s_ab = values[1], values[2], float(values[4])
            if protein_b < protein_a:
                raise Exception('protein pairs should be in alphabetical order')
            proteins.add(protein_a)
            proteins.add(protein_b)
            s_values[(protein_a, protein_b)] = s_ab
    proteins = sorted(list(proteins))

    # Assert there are n choose 2 combinations, where n = len(proteins)
    assert len(s_values) == (len(proteins) - 1) * len(proteins) / 2

    column_labels = row_labels = proteins
    data = np.zeros((len(proteins), len(proteins)))
    for i, protein_a in enumerate(proteins):
        for j, protein_b in enumerate(proteins):
            if protein_a == protein_b:
                data[i,j] = 0
            else:
                data[i,j] = s_values[tuple(sorted((protein_a, protein_b)))]
    fig, ax = plt.subplots()
    heatmap = ax.pcolor(data, cmap=plt.cm.Blues)

    fig.suptitle(title, fontsize=13)

    ax.set_xticks(np.arange(data.shape[0])+0.5, minor=False)
    ax.set_yticks(np.arange(data.shape[1])+0.5, minor=False)

    ax.invert_yaxis()
    ax.xaxis.tick_top()

    ax.set_xticklabels(row_labels, minor=False)
    ax.set_yticklabels(column_labels, minor=False)

    plt.savefig(output)

if __name__ == '__main__':
    main()
