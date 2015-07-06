# Converts the original flu genes file (original_flu_genes.csv) from Uniprot
# Accession IDs to Gene IDs, and only retains human genes.
#
# Usage: python scripts/01_filter_flu_genes.py > data/flu_gene_ids.tsv

from convert_uniprot import gene_ids_from_uniprot_accessions

print '# Flu Proteins/Genes'
print '# PROTEIN  GENE'
with open('data/original_flu_genes.csv') as f:
    proteins = []
    uniprots = []
    for line in f:
        values = line.split(',')
        entry_name = values[4]
        if entry_name.strip().endswith('_HUMAN'):
            proteins.append(values[0])
            uniprots.append(values[1])

    gene_ids = gene_ids_from_uniprot_accessions(uniprots)

    for protein, gene_id in zip(proteins, gene_ids):
        if gene_id:
            print (protein + '\t' + ';'.join(gene_id))
