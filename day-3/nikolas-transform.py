#!/usr/bin/env python

import sys
import csv
from ete3 import Tree
from collections import Counter


def loadNames():
    names = {}

    with open("../niklas_problem/annotations/all.csv") as fp: # Why do it like this
        reader = csv.reader(fp, delimiter=",")

        for fields in reader:
            oldName, newName = fields[8], fields[6]

            if oldName == "Protein product":
                assert fields[0] == "#Name"
            else:
                names[oldName] = newName

    return names


def relabelTree(filename, newNames):

    missCount = 0
    speciesCounts = Counter()



    tree = Tree(filename, format=1)
    for leaf in tree.get_leaves():
        fields = leaf.name.split("_")
        assert len(fields) == 4
        species = fields[0] + "_" + fields[1]
        protein = fields[2] + "_" + fields[3]

        if protein in newNames:
            leaf.name = species + "_" + newNames[protein]
        else:
            missCount += 1

        speciesCounts[species] += 1

    newFilename = filename + ".new"
    tree.write(format=1, outfile=newFilename)

    if missCount > 0:
        print(f"Hey, there were {missCount} misses in {filename}!")

    return speciesCounts

# test
def main():
    newNames = loadNames()
    speciesCounts = Counter()

    for filename in sys.argv[1:]:
        # print("Processing", filename)
        speciesCounts += relabelTree(filename, newNames)

    print("The species counts are:")
    for species, count in speciesCounts.items():
        print(f"   {species}: {count}")


if __name__ == '__main__':
    main()
