from xml.dom import minidom
import requests

def convert_to_gene_symbols(input, output):
    gene_ids = set()
    with open(input) as f:
        for line in f.read().splitlines():
            if line.startswith('#'):
                continue
            gene_ids.update(line.split('\t')[1].split(';'))
    gene_ids = list(gene_ids)

    gene_symbols = {}
    groups = [gene_ids[i:i+200] for i in range(0, len(gene_ids), 200)]
    for group in groups:
        url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id=%s' % ','.join(group)
        data = requests.get(url).content
        xmldoc = minidom.parseString(data)
        itemlist = xmldoc.getElementsByTagName('DocumentSummary')
        for item in itemlist:
            gene_id = item.attributes['uid'].value
            gene_symbol = item.getElementsByTagName('Name')[0].childNodes[0].nodeValue
            gene_symbols[gene_id] = gene_symbol

    with open(output, 'wb') as out:
        for gene_id in gene_ids:
            gene_symbol = gene_symbols[gene_id]
            out.write('\t'.join((gene_id, gene_symbol)))
            out.write('\n')

convert_to_gene_symbols('data/flu_gene_ids.tsv', 'data/flu_gene_symbols.tsv')
convert_to_gene_symbols('data/hek_gene_ids.tsv', 'data/hek_gene_symbols.tsv')
convert_to_gene_symbols('data/jurkat_gene_ids.tsv', 'data/jurkat_gene_symbols.tsv')
