from collections import defaultdict
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import math
import numpy as np
import pandas as pd
import seaborn as sns

def generate_stripplot(max_distance, category_filter, title, output):
    protein_pairs = ('HA-GP160', 'M2-GP160', 'M2-VPR', 'NA-GP160', 'NA-VPU')
    ontologies = ('biological_process', 'cellular_component', 'molecular_function')

    labels = {}
    p_values = defaultdict(list)
    counts = defaultdict(int)

    for protein_pair in protein_pairs:
        for ontology in ontologies:
            file = 'output/ontology/' + protein_pair + '_' + ontology + '_up-to-' + str(max_distance) + '.tsv'
            with open(file) as f:
                for line in f.read().splitlines()[2:]:
                    label, ref, ID, url, p_value = line.split('\t')
                    p_value = float(p_value)
                    if p_value > 0.05 or p_value == 0:
                        continue
                    labels[ID] = label
                    p_values[ID].append(float(p_value))
                    counts[ID] += 1

    p_value_averages = p_values
    for ID in p_values:
        p_value_averages[ID] = np.mean(p_values[ID])
    IDs = sorted(sorted(labels.keys()), key=lambda x: category_filter(x, labels))

    y = np.array(map(lambda x: counts[x], IDs))
    x = np.array(map(lambda x: -math.log(p_value_averages[x], 10), IDs))
    colors = map(lambda x: 'r' if category_filter(x, labels) else 'b', IDs)
    areas = map(lambda x: 80 if category_filter(x, labels) else 20, IDs)

    for ID in sorted(sorted(filter(lambda x: category_filter(x, labels), IDs), key=p_value_averages.get), key=counts.get, reverse=True):
        label, count, p_value_avg = labels[ID], counts[ID], p_value_averages[ID]
        print '\t'.join((ID, label, str(counts[ID]), str(p_value_avg)))
    print

    plt.clf()
    plt.scatter(x, y, c=colors, s=areas)
    plt.gca().set_xlabel("-log(p-value)")
    plt.gca().set_ylabel("count")
    plt.suptitle(title)
    plt.savefig(output)

def is_ER(ID, labels):
    label = labels[ID]
    return ('ER' in label or 'endoplasmic reticulum' in label.lower())

def is_localization(ID, labels):
    label = labels[ID]
    return 'protein localization' in label.lower()

def is_localization_or_ER(ID, labels):
    return is_ER(ID, labels) and is_localization(ID, labels)

generate_stripplot(1, is_ER, 'Endoplasmic Reticulum p-Values (d <= 1)',
                   'output/figures/ER_stripplot_up-to-1.png')
generate_stripplot(1, is_localization, 'Protein Localization p-Values (d <= 1)',
                   'output/figures/protein-localization_stripplot_up-to-1.png')
generate_stripplot(1, is_localization_or_ER, 'Endoplasmic Reticulum & Protein Localization p-Values (d <= 1)',
                   'output/figures/protein-localization_ER_up-to-1.png')
generate_stripplot(0, is_ER, 'Endoplasmic Reticulum p-Values (d = 0)',
                   'output/figures/ER_stripplot_up-to-0.png')
generate_stripplot(0, is_localization, 'Protein Localization p-Values (d = 0)',
                   'output/figures/protein-localization_stripplot_up-to-0.png')
generate_stripplot(0, is_localization_or_ER, 'Endoplasmic Reticulum & Protein Localization p-Values (d = 0)',
                   'output/figures/protein-localization_ER_up-to-0.png')
