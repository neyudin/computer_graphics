from sys import argv
from statistics import mean

values = []
with open(argv[1]) as f:
    flag = True
    for line in f:
        if flag:
            flag = False
            continue
        values.append(float(line.split()[-1]))

print('Mean precision: {}'.format(mean(values)))
