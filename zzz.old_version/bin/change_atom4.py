#!/usr/bin/env pythonw
from visual import *
from numpy import *
import lift
import sys

cell = lift.Cell()
if		len(sys.argv) == 1:	cell.load('input.lat', 'c')
elif	len(sys.argv) == 2:	cell.load(sys.argv[1], 'o')
elif	len(sys.argv) == 3:	cell.load(sys.argv[1], sys.argv[2])


for i in range(len(cell.superlattice)):
    AtomIndex, UCIndex, Name, pos = cell.superlattice[i]

    x = pos[0,0]
    y = pos[0,1]
    z = pos[0,2]
    #if x > 1.9 and x< 7.1 and y>1.9 and y<7.1:

    if x > 3.6 and x <5.6 and y > 3.6 and y <5.6 :
        Name = Name.replace("Fe", "Mn")
        Name = Name.replace("Bi", "LS")

        #Name = Name.replace("Mn", "Fe")
        #Name = Name.replace("LS", "Bi")

    #if x < 0.6:
    #    Name = Name.replace("Fe", "VA")
    #    Name = Name.replace("Bi", "VA")
    #    Name = Name.replace("O",  "VA")

    cell.superlattice[i] = AtomIndex, UCIndex, Name, pos


cell.save()
