import re
import requests

def gene_ids_from_uniprot_accessions(uniprot_ids):
    url = 'http://biodbnet.abcc.ncifcrf.gov/db/db2dbRes.php'

    # Split Uniprot IDs into groups of 300
    groups = [uniprot_ids[i:i+300] for i in range(0, len(uniprot_ids), 300)]
    results = []

    for group in groups:
        payload = {
            'input': 'Uniprot Accession',
            'outputs[]': 'Gene ID',
            'taxonId': 'optional',
            'idList': '\n'.join(group),
            'hasComma': 'no',
            'removeDupValues': 'no',
            'expandTaxId': 'yes',
            'request': 'db2db'
        }
        resp = requests.post(url, data=payload).content
        link = re.search(r"href='(\/db\/dbResFileTxt\.php.*)'", resp).group(1)
        link = 'http://biodbnet.abcc.ncifcrf.gov' + link

        conversion = requests.get(link).content
        results.extend(map(_gene_ids_from_row, conversion.splitlines()[1:]))

    assert len(uniprot_ids) == len(results)
    return results

def _gene_ids_from_row(row):
    gene_ids = row.split('\t')[1]
    if gene_ids == '-':
        return None
    return tuple(gene_ids.split('; '))
