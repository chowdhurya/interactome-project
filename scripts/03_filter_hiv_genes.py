# Converts HIV genes from Uniprot Accession to Gene ID
# Turns:
#     - original_hek_genes.csv into hek_gene_ids.tsv
#     - original_jurkat_genes.csv into jurkat_gene_ids.tsv
# Usage: python scripts/03_filter_hiv_genes.py

from convert_uniprot import gene_ids_from_uniprot_accessions

def convert_file(input, output):
    proteins = []
    uniprots = []
    with open(input) as f:
        for line in f.read().splitlines():
            values = line.split(',')
            proteins.append(values[0])
            uniprots.append(values[1])

    gene_ids = gene_ids_from_uniprot_accessions(uniprots)

    with open(output, 'wb') as out:
        out.write('#PROTEIN\tGENE ID\n')
        for protein, gene_id in zip(proteins, gene_ids):
            if gene_id:
                out.write(protein + '\t' + ';'.join(gene_id) + '\n')

convert_file('data/original_jurkat_genes.csv', 'data/jurkat_gene_ids.tsv')
convert_file('data/original_hek_genes.csv', 'data/hek_gene_ids.tsv')
