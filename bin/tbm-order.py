#!/usr/bin/env python
#-------------------------------------------------------------|
#| Copyright (C) 2016 Yuan-Yen Tai, Hongchul Choi,            |
#|                    Jian-Xin Zhu                            |
#|                                                            |
#| This file is distributed under the terms of the BSD        |
#| Berkeley Software Distribution. See the file `LICENSE' in  |
#| the root directory of the present distribution, or         |
#| https://en.wikipedia.org/wiki/BSD_licenses.                |
#|                                                            |
#|-------------------------------------------------------------
import os
import sys
import scipy
import numpy as np

from scipy import spatial

def MatToStr(valMat):
	valStr = ""
	for val in np.array(valMat)[0]:
		valStr += str(val)+","
	return valStr[0:-1]




class AtomOrder:
	def __init__(self, atomStringList, _atomInfo):
		#----------------------------------
		# Parse the following structure.
		# >>> 1   7    Fe   [[  0.5         0.5         0.5       ]]
		#----------------------------------
		self.atomInfo = _atomInfo
		self.atomIndex = int(atomStringList[1])
		self.atomName = atomStringList[3]
		self.atomPos = map(float, atomStringList[5:8])

		self.orderMap = {}
		self.orderSequence = []

	def appendOrder(self, orderString):
		if len(orderString) != 2:
			return
		name = orderString[0]
		name = name.replace(" ", "")
		order = scipy.mat(map(float, orderString[1].split(',')))
		self.orderMap[name] = order
		self.orderSequence.append(name)

	def getOrderString(self):
		orderString = ""
		for orderName in self.orderSequence:
			orderString += " "+orderName + " = " +  MatToStr(self.orderMap[orderName]) + "\n"

		return orderString





class   LatticeOrder:
	def __init__(self, _filename):
		self.filename = _filename
		f = open(self.filename)

		self.atomOrderList = []
		self.atomPos = []
		for line in f.readlines():
			spline = line.split()

			if spline[0] == '>>>':
				self.atomOrderList.append(AtomOrder(spline, line))
				self.atomPos.append((
					self.atomOrderList[-1].atomPos[0],
					self.atomOrderList[-1].atomPos[1],
					self.atomOrderList[-1].atomPos[2]
				))
			else:
				order = line.split('=')
				self.atomOrderList[-1].appendOrder(order)

		f.close()
		# Using the atom position to build the rtree.
		self.atomPos = np.array(self.atomPos)
		self.tree3D = spatial.KDTree( self.atomPos )

	def getOrder(self, key, pos):
		# ---------------------------------------------------------------
		# To get the correspond order for an atom of a given position
		# ---------------------------------------------------------------
		atomName, orderName = key
		result = self.tree3D.query(np.array([pos]))

		if result[0] < 0.01 :
			atom = self.atomOrderList[int(result[1])]
			atomOrder = self.atomOrderList[int(result[1])].orderMap

			if atomName == "":
				atomName = atom.atomName

			if (orderName in atomOrder) and atomName == atom.atomName:
				return True,atomOrder[orderName].copy()

		return False,np.mat([0.0])

	def setOrder(self, key, pos, order):
		# ---------------------------------------------------------------
		# To set the correspond order for an atom of a given position
		# ---------------------------------------------------------------
		atomName, orderName = key
		result = self.tree3D.query(np.array([pos]))

		if result[0] < 0.01 :
			atom = self.atomOrderList[int(result[1])]
			atomOrder = self.atomOrderList[int(result[1])].orderMap

			if atomName == "":
				atomName = atom.atomName

			if (orderName in atomOrder) and atomName == atom.atomName:
				print atom.atomInfo,
				print orderName, " = ", MatToStr( np.mat( order ) )
				self.atomOrderList[int(result[1])].orderMap[orderName] = np.mat(order)

		#if result[0] < 0.01 and orderName in self.atomOrderList[int(result[1])].orderMap:
		#	atom = self.atomOrderList[int(result[1])]
		#	atomOrder = self.atomOrderList[int(result[1])]

		#	if atomName == "":
				atomName = atom.atomName


	def save(self, _filename = ""):
		fname = ""
		fname == self.filename
		if _filename != "":
			fname = _filename+".ord"

		f = open(fname, 'w')
		for atom in self.atomOrderList:
			f.write(atom.atomInfo)
			f.write(atom.getOrderString())
		f.close()





if __name__ == "__main__":

	inputLatticeFileName = "sys.argv[1]"

	print sys.argv

	content = """#!/usr/bin/env python
#|---------------------------------------------------------------
#| Modify this file to manipulate the orders in following ways:
#| 1. Access the order:
#|         LatOrder.getOrder( key, pos)
#| 2. Set new value to the order:
#|         LatOrder.setOrder( key, pos, order)
#| 3. Save the new order:
#|         LatOrder.save( filename = "" )
#|---------------------------------------------------------------
import imp
import sys
import numpy as np

### Loading LatticeOrder class from the TBM3 bin/ path.
foo = imp.load_source('LO', '"""+sys.argv[0]+"""')
LatOrder = foo.LatticeOrder("""+inputLatticeFileName+"""+'.ord')

#### Query the correspond order in the given atom position.
#found, order = LatOrder.getOrder(key=('Fe','@:cspin'), pos=[0, 0, 0])

#### Setup the correspond order in the given atom position.
#LatOrder.setOrder(key=('Fe','@:cspin'), pos=[0.5, 0.5, 0.5,], order=[1,2,3])

#### Save the order in the original input file.
#LatOrder.save()
"""
	filename = "orderAnalyzer.py"
	if len(sys.argv) == 2:
		filename = sys.argv[1]+".py"

	f = open(filename,'w')
	f.write(content)
	f.close()

	os.system("chmod +x "+filename)

