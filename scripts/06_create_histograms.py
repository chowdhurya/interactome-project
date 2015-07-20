import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.pyplot import savefig
from collections import defaultdict, Counter
from itertools import groupby
import numpy as np
import pandas as pd
import re
import seaborn as sns
import StringIO

def parse_protein(protein):
    m = re.search(r'(?P<protein_name>.+) \[(?P<disease>.*)\]', protein)
    return {'protein': m.group('protein_name'), 'disease': m.group('disease')}

print """<!DOCTYPE html>
<html lang="en">
  <head>
    <meta charset="utf-8">
    <title>Flu/Jurkat Gene Distance Histograms</title>
    <link href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.5/css/bootstrap.min.css" rel="stylesheet">
  </head>
  <body>"""

with open('output/flu-jurkat_distances.txt') as f:
    data = [list(group) for k, group in groupby(f.read().splitlines(), bool) if k][1:]
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
            distances.append(int(re.search(r'\d+ -> \d+ = (\d+)', dist_data).group(1)))

        plt.clf()
        sns.countplot(data=pd.DataFrame(distances, columns=['distance']), x='distance')
        plt.suptitle(protein_a['protein'] + ' [Flu] / ' + protein_b['protein'] + ' [Jerkat] Gene-Gene Distances\ns_AB = ' + str(s_AB)   , fontsize=13)
        sio = StringIO.StringIO()
        plt.savefig(sio, format='png')

        print '  <div class="row">'
        print '    <div class="col-xs-6 col-xs-offset-3">'
        print '      <h2 class="text-center">%s [Flu] and %s [Jurkat]</h2>' % (protein_a['protein'], protein_b['protein'])
        print '      <img src="data:image/png;base64,%s"/>' % sio.getvalue().encode("base64").strip().replace('\n', '')
        print '    </div>'
        print '  </div>'

print """  </body>
</html>"""
