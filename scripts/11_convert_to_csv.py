import os

for file in os.listdir('output/ontology'):
    if not file.endswith('.tsv'):
        continue
    with open('output/ontology/' + file) as f:
        with open('output/ontology/csv/' + file[:-4] + '.csv', 'wb') as out:
            for line in f.read().splitlines():
                out.write(','.join(map(lambda x: '"%s"' % x, line.split('\t'))))
                out.write('\n')
