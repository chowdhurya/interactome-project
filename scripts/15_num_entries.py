import requests

protein_pairs = ('HA-GP160', 'M2-GP160', 'M2-VPR', 'NA-GP160', 'NA-VPU')
ontologies = ('biological_process', 'cellular_component', 'molecular_function')
max_distance = 1

num_entries_cache = {}
def _get_num_entries(ID):
    if ID in num_entries_cache:
        return
    url = 'http://geneontology-golr.stanford.edu/solr/select?&wt=json&fq=regulates_closure:"%s"&fq=document_category:"annotation"&q=*:*' % (ID,)
    num_entries = requests.get(url).json()['response']['numFound']
    num_entries_cache[ID] = num_entries

for protein_pair in protein_pairs:
    for ontology in ontologies:
        file = 'output/ontology/' + protein_pair + '_' + ontology + '_up-to-' + str(max_distance) + '.tsv'
        with open(file) as f:
            for i, line in enumerate(f.read().splitlines()[2:]):
                print i
                label, ref, ID, url, p_value = line.split('\t')
                _get_num_entries(ID)

with open('data/num_entries.tsv', 'wb') as out:
    out.write('ID\tNum Entries\n')
    for ID in sorted(num_entries_cache.keys()):
        num_entries = num_entries_cache[ID]
        if num_entries == 0:
            continue
        out.write('\t'.join((ID, str(num_entries))))
        out.write('\n')
