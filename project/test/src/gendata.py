import os
import shutil
from sys import argv
from random import random

if os.path.isdir('./data/multiclass'):
    shutil.rmtree('./data/multiclass')

os.mkdir('./data/multiclass')
os.mkdir('./data/multiclass/train')
os.mkdir('./data/multiclass/test')

data = {}
dirs = {}
# load data
with open('./data/whole_labels.txt') as base:
    for line in base:
        filename, label = line.strip().split()
        filename = filename[filename.find('/')+1:]
        dirs[label] = filename[:filename.find('/')]
        try:
            data[label].append(filename)
        except KeyError:
            data[label] = [filename]

ratio = float(argv[1])

with open('./data/multiclass/train_labels.txt', 'w') as train:
    with open('./data/multiclass/test_labels.txt', 'w') as test:
        for label, filenames in data.items():

            os.mkdir('./data/multiclass/train/{}'.format(dirs[label]))
            os.mkdir('./data/multiclass/test/{}'.format(dirs[label]))

            for i, filename in enumerate(filenames):
                if random() < ratio:
                    train.write('train/{} {}\n'.format(filename, label))
                    shutil.copy('../template/data/multiclass/train/{}'.format(filename), './data/multiclass/train/{}'.format(filename))
                else:
                    test.write('./test/{} {}\n'.format(filename, label))
                    shutil.copy('../template/data/multiclass/train/{}'.format(filename), './data/multiclass/test/{}'.format(filename))

