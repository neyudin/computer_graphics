#!/usr/bin/python

import sys
import os

def load_labels(filename):
    lines = open(filename, 'r').readlines()
    dict_labels = {}
    labels = []
    for line in lines:
        parts = line.split(' ')
        if not dict_labels.__contains__(parts[0]):
            dict_labels[parts[0]] = parts[1]
            labels.append(parts[1])
    return labels


def test_labels(gt_labels, predicted_labels):
    if (len(gt_labels) != len(predicted_labels)):
        print "Error! Files with predicted and ground truth labels " \
              "have different number of samples."
        return
    if (len(gt_labels) == 0):
        print "Error! Dataset is empty."
        return

    correct_len = 0
    for gt_elem, pred_elem in zip(gt_labels, predicted_labels):
        if gt_elem == pred_elem:
            correct_len += 1
    precision = float(correct_len) / len(gt_labels)
    sys.stdout = open(sys.argv[3], 'a')
    print "Precision: %f" % precision

if len(sys.argv) != 4:
    print 'Usage: %s <ground_truth.txt> <program_output.txt> <stdout.txt>' % sys.argv[0]
    sys.exit(0)

gt_labels = load_labels(sys.argv[1])
predicted_labels = load_labels(sys.argv[2])

test_labels(gt_labels, predicted_labels)
