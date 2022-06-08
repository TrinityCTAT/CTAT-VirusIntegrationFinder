#!/usr/bin/env python

import sys, os, re


clusters = list()

with open("virus_db.fasta.cdhit.clstr") as fh:

    current_cluster = list()
    for line in fh:
        line = line.rstrip()
        if line[0] == ">":
            current_cluster = list()
            clusters.append(current_cluster)
        else:
            vals = re.split("\\s+", line)
            acc = vals[2]
            acc = acc[1:-3]
            current_cluster.append(acc)

for cluster in clusters:
    rep_acc = cluster[0]
    if len(cluster) > 1:
        hpv_accs = list(filter(lambda x: re.match('HPV', x), cluster))
        if hpv_accs:
            rep_acc = hpv_accs[0]

    print(rep_acc)

sys.exit(0)





