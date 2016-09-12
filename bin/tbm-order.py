#!/usr/bin/env python

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

	def getOrder(self, orderName, pos):
		# ---------------------------------------------------------------
		# To get the correspond order for an atom of a given position
		# ---------------------------------------------------------------
		result = self.tree3D.query(np.array([pos]))
		if result[0] < 0.01 :
			order = self.atomOrderList[int(result[1])].orderMap
			if orderName in order:
				return True,order[orderName].copy()

		return False,np.mat([0.0])

	def setOrder(self, orderName, pos, order):
		# ---------------------------------------------------------------
		# To set the correspond order for an atom of a given position
		# ---------------------------------------------------------------
		result = self.tree3D.query(np.array([pos]))
		if result[0] < 0.01 and orderName in self.atomOrderList[int(result[1])].orderMap:
			atomOrder = self.atomOrderList[int(result[1])]
			print atomOrder.atomInfo
			self.atomOrderList[int(result[1])].orderMap[orderName] = np.mat(order)

	def save(self):
		f = open(self.filename, 'w')
		for atom in self.atomOrderList:
			f.write(atom.atomInfo)
			f.write(atom.getOrderString())
		f.close()





if __name__ == "__main__":
	content = """#!/usr/bin/env python
import imp
import sys
import numpy as np

foo = imp.load_source('LO', '"""+sys.argv[0]+"""')
LatOrder = foo.LatticeOrder(sys.argv[1]+".ord")
#found, order = LatOrder.getOrder('@:cspin', pos=[0, 0, 0])
#LatOrder.setOrder('@:cspin', pos=[0.5, 0.5, 0.5,] order=[1,2,3])
#LatOrder.save()
"""
	filename = "orderAnalyzer.py"
	if len(sys.argv) == 2:
		filename = sys.argv[1]+".py"

	f = open(filename,'w')
	f.write(content)
	f.close()

	os.system("chmod +x "+filename)

