#!/usr/bin/env python

import numpy as np
from math import *
import bisect

# The basic class data structure for the coord_search
class coord3D:
	radius = 0.001
	def __init__(self, pos, data):
		self.x = pos[0,0]
		self.y = pos[0,1]
		self.z = pos[0,2]
		self.data = data

	def __lt__(self, other):
		if abs(self.x-other.x) < coord3D.radius:
			if abs(self.y-other.y) < coord3D.radius:
				if abs(self.z-other.z) < coord3D.radius: return False
				else : return self.z < other.z
			else: return self.y < other.y
		else: return self.x < other.x

	def pos(self):
		return np.mat([[self.x,self.y,self.z]])

# The function use to calculate the distance between two points.
def distance(p0, p1):
	x = abs(p0.x - p1.x)
	y = abs(p0.y - p1.y)
	z = abs(p0.z - p1.z)
	return sqrt(x**2+y**2+z**2)

class coord_search:
	def __init__(self, radius = 0.001):
		coord3D.radius = radius
		self.list = []

	def add_point(self, pos, data):
			bisect.insort(self.list, coord3D(pos, data))

	def Print(self):
		for i in self.list:
			print i.x, i.y, i.z

	def search(self, point):
		index_found = bisect.bisect_left(self.list, coord3D(point,0))

		def if_return_index(index, pp):
			p0 = self.list[index_found]
			p1 = coord3D(pp,0)
			if distance(p0,p1) < coord3D.radius : return index, p0
			else: return -1, p0

		return if_return_index(index_found, point)

if __name__ == '__main__':

	A = coord_search()
	A.add_point(np.mat([-1,2,3]), 1)
	A.add_point(np.mat([1,3,3]), 2)
	A.add_point(np.mat([1,4,3]), 3)
	A.add_point(np.mat([2,2,3]), 4)
	A.Print()

	ii, at=A.search(np.mat([0.99999,3,3]))
	print ii, at.pos()



