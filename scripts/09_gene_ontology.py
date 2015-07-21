from itertools import groupby
import re
import requests

# Read in gene symbols
flu_gene_symbols = {}
jurkat_gene_symbols = {}
with open('data/flu_gene_symbols.tsv') as f:
    flu_gene_symbols.update(map(lambda x: x.split('\t'), f.read().splitlines()))
with open('data/jurkat_gene_symbols.tsv') as f:
    jurkat_gene_symbols.update(map(lambda x: x.split('\t'), f.read().splitlines()))

def gene_id_to_symbol(id):
    if id in flu_gene_symbols:
        return flu_gene_symbols[id]
    if id in jurkat_gene_symbols:
        return jurkat_gene_symbols[id]
    raise Exception('ID not found')

def parse_protein(protein):
    m = re.search(r'(?P<protein_name>.+) \[(?P<disease>.*)\]', protein)
    return {'protein': m.group('protein_name'), 'disease': m.group('disease')}

with open('output/flu-jurkat_distances.txt') as f:
    data = [list(group) for k, group in groupby(f.read().splitlines(), bool) if k][1:]

    i = 0
    for data_set in data:
        m = re.search(r'\d. Proteins: (?P<protein_a>.*) and (?P<protein_b>.*)', data_set[0])
        protein_a, protein_b = parse_protein(m.group('protein_a')), parse_protein(m.group('protein_b'))
        if protein_a['disease'] == protein_b['disease']:
            continue
        if protein_a['disease'] == 'Jurkat':
            protein_a, protein_b = protein_b, protein_a
        assert protein_a['disease'] == 'Flu' and protein_b['disease'] == 'Jurkat'

        distances = []
        for dist_data in data_set[3:]:
            distances.append(re.search(r'(\d+) -> (\d+) = (\d+)', dist_data).groups())
        distances = map(lambda x: (gene_id_to_symbol(x[0]), gene_id_to_symbol(x[1]), int(x[2])), distances)

    	for max_distance in (0, 1):
    		for ontology in ('biological_process', 'molecular_function', 'cellular_component'):
    			payload = {
    				'ontology': ontology,
    				'species': 'HUMAN',
    				'correction': 'bonferroni',
    				'format': 'json',
    				'input': '\n'.join(map(lambda x: x[0], filter(lambda x: x[2] <= max_distance, distances)))
    			}
    			results = requests.get('http://pantherdb.org/webservices/go/overrep.jsp', params=payload).json()['results']['result']
    			results = filter(lambda x: x['pValue'] < 0.5, results)
    			filename = 'output/ontology/%s-%s_%s_up-to-%d.tsv' % (protein_a['protein'], protein_b['protein'], ontology, max_distance)
    			with open(filename, 'wb') as out:
    				out.write('%s (Flu) / %s (Jurkat) %s Ontology, dist <= %d\n' % (protein_a['protein'], protein_b['protein'], ontology, max_distance))
    				out.write('\t'.join(('Label', 'REF', 'ID', 'URL', 'p-value')))
    				out.write('\n')
    				for result in results:
    					label = result['term']['label']
    					ref = result['number_in_reference']
    					go_id = result['term'].get('id', '')
    					url = ('http://amigo.geneontology.org/amigo/term/%s' % go_id) if go_id else ''
    					p_value = result['pValue']
    					out.write('\t'.join((label, str(ref), go_id, url, str(p_value))))
    					out.write('\n')
    			
    	# Do it for first 5 protein pairs
    	i += 1
    	if i == 5:
    		break