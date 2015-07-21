from itertools import groupby
import re

# Read in gene symbols
flu_gene_symbols = {}
jurkat_gene_symbols = {}
with open('data/flu_gene_symbols.tsv') as f:
    flu_gene_symbols.update(map(lambda x: x.split('\t'), f.read().splitlines()))
with open('data/jurkat_gene_symbols.tsv') as f:
    jurkat_gene_symbols.update(map(lambda x: x.split('\t'), f.read().splitlines()))

def gene_id_to_symbol(id):
    if id in flu_gene_symbols and id in jurkat_gene_symbols:
        return '%s (%s) [Both]' % (flu_gene_symbols[id], id)
    if id in flu_gene_symbols:
        return '%s (%s) [Flu]' % (flu_gene_symbols[id], id)
    if id in jurkat_gene_symbols:
        return '%s (%s) [Jurkat]' % (jurkat_gene_symbols[id], id)
    raise Exception('ID not found')

def parse_protein(protein):
    m = re.search(r'(?P<protein_name>.+) \[(?P<disease>.*)\]', protein)
    return {'protein': m.group('protein_name'), 'disease': m.group('disease')}

with open('output/flu-jurkat_distances.txt') as f:
    with open('output/flu-jurkat_crossed_symbol_distances.txt', 'wb') as out:
        data = [list(group) for k, group in groupby(f.read().splitlines(), bool) if k][1:]
        i = 1
        for data_set in data:
            m = re.search(r'\d. Proteins: (?P<protein_a>.*) and (?P<protein_b>.*)', data_set[0])
            protein_a, protein_b = parse_protein(m.group('protein_a')), parse_protein(m.group('protein_b'))
            if protein_a['disease'] == protein_b['disease']:
                continue
            if protein_a['disease'] == 'Jurkat':
                protein_a, protein_b = protein_b, protein_a
            assert protein_a['disease'] == 'Flu' and protein_b['disease'] == 'Jurkat'
            s_AB = float(re.search(r's_AB = (?P<s_AB>.*)', data_set[2]).group(1))

            distances = []
            for dist_data in data_set[3:]:
                distances.append(re.search(r'(\d+) -> (\d+) = (\d+)', dist_data).groups())
            distances = map(lambda x: (gene_id_to_symbol(x[0]), gene_id_to_symbol(x[1]), x[2]), distances)

            out.write('%d. Proteins: %s [Flu] and %s [Jurkat] \n' % (i, protein_a['protein'], protein_b['protein']))
            out.write('s_AB = %f\n' % s_AB)
            for distance in distances:
                out.write('%s -> %s = %s\n' % distance)
            out.write('\n')
            i += 1
