#!/usr/bin/env python
import sys
import scipy
import numpy as np
from scipy import spatial


if len(sys.argv) == 2:
    fileStorage = []
    f = open(sys.argv[1], 'r')

    header = ""
    atomList = []
    posDic = {}
    for line in f.readlines():
        line = line.replace('\t', ' ')
        parser = line.split()
        if len(parser) > 0:
            if parser[0][0] == "#":
                header = parser[0]


        if header != "#Atoms":
            fileStorage.append(line)
        else:
            atom = []
            if len(parser) == 4:
                for elem in parser:
                    atom.append(float(elem))

                atomList.append(atom)

    print atomList
