#!/usr/bin/env python

import sys, csv

should = False
dircnt = {}

def read_csv(dst, path):
    with open(path, 'r') as f:
        r = csv.reader(f, delimiter='\t')
        last0 = ""
        cnt = 0
        for row in r:
            if row[0] == last0:
                if not should:
                    dircnt[row[0]] += 1
                cnt += 1
            else:
                if not should:
                    dircnt[row[0]] = 1
                cnt = 1
            granica = int(sys.argv[3])
            if should:
                d = 0
                if row[0] in dircnt:
                    d = dircnt[row[0]]
                granica = min(granica, d)
            if (cnt <= granica):
                dst.append((row[0], row[1]))
            last0 = row[0]

a = []
read_csv(a, sys.argv[1])

b = []
should = True
read_csv(b, sys.argv[2])
print ("testing first " + str(sys.argv[3]))
print ("blast out first")
print(len(set(a)))
print(len(set(b)))
print(len(set(a) & set(b)))

print(float(len(set(a) & set(b))) / len(set(a) | set(b)))
