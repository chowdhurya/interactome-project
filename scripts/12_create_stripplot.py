from collections import defaultdict
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

protein_pairs = ('HA-GP160', 'M2-GP160', 'M2-VPR', 'NA-GP160', 'NA-VPU')
ontologies = ('biological_process', 'cellular_component', 'molecular_function')

labels = {}
p_values = defaultdict(list)
counts = defaultdict(int)

for protein_pair in protein_pairs:
    for ontology in ontologies:
        file = 'output/ontology/' + protein_pair + '_' + ontology + '_up-to-1.tsv'
        with open(file) as f:
            for line in f.read().splitlines()[2:]:
                label, ref, ID, url, p_value = line.split('\t')
                p_value = float(p_value)
                if p_value > 0.05:
                    continue
                labels[ID] = label
                p_values[ID].append(float(p_value))
                counts[ID] += 1

p_value_averages = p_values
for ID in p_values:
    p_value_averages[ID] = np.mean(p_values[ID])

IDs = sorted(labels.keys())
data = pd.DataFrame.from_records(map(lambda ID: (p_value_averages[ID], 'x' + str(counts[ID]).zfill(2)), IDs),
                                 columns=('p_value_avg', 'count'))

sns.stripplot(x='p_value_avg', y='count', data=data.sort('count', ascending=False))
plt.savefig('stripplot.png')
