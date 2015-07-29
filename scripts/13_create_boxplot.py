from collections import defaultdict
import math
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
                if p_value > 0.05 or p_value == 0:
                    continue
                labels[ID] = label
                p_values[ID].append(float(p_value))
                counts[ID] += 1

def is_ER(ID):
    label = labels[ID]
    return 'ER' in label or 'endoplasmic reticulum' in label.lower()

p_value_averages = p_values
for ID in p_values:
    p_value_averages[ID] = np.mean(p_values[ID])
IDs = sorted(sorted(labels.keys()), key=is_ER)

count_strings = sorted(set(map(lambda ID: 'x' + str(counts[ID]).zfill(2), IDs)), reverse=True)
data = pd.DataFrame.from_records(map(lambda ID: ('x' + str(counts[ID]).zfill(2), math.log(p_value_averages[ID])), IDs),
                              columns=('count', 'ln(average p-value)'))

sns.boxplot(x='ln(average p-value)', y='count', data=data, order=count_strings)
plt.savefig('output/figures/boxplot.png')
