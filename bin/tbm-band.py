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
from numpy import *
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid

# The recommended way to use wx with mpl is with the WXAgg backend. 
import matplotlib


class Band:
	def __init__(self,filename,shift=0):

		if filename[-3:] != "ban":
			filename = filename+".ban"

		f = open(filename,"r")

		self.pointLabel = []
		self.arrayKpoin = []
		self.kx = []
		self.ky = []
		self.kz = []
		self.Elist= []

		for line in f.readlines():
			sp = line.split()
			self.pointLabel.append(sp[0])
			kx = float(sp[2])
			ky = float(sp[3])
			kz = float(sp[4])
			self.arrayKpoin.append(array([kx,ky,kz]))
			En = []
			for ss in sp[7:-1] :
				En.append(float(ss)-shift)

			self.Elist.append(En)

		self.Elist = mat(self.Elist)

	def plot(self):
		Ek = self.Elist.T.tolist()

		Klength = 0
		Kpath = []
		Kpath.append(0)

		KVerticalLinePos = []

		x_label = []

		plt.axhline(y=0, linestyle=":", linewidth=1, color='k')

		for i in xrange(1,len(self.arrayKpoin)):
			kp_0 = self.arrayKpoin[i-1]
			kp_1 = self.arrayKpoin[i]
			dkp = kp_1 - kp_0
			ds = sqrt(dkp[0]**2 + dkp[1]**2 + dkp[2]**2)
			Klength += ds
			Kpath.append(Klength)

			if self.pointLabel[i-1] != '-':
				KVerticalLinePos.append(Kpath[i-1])
				x_label.append(self.pointLabel[i-1])
			else:
				x_label.append(" ")


			if i == len(self.arrayKpoin)-1 :
				x_label.append(self.pointLabel[i])

		x_ticks = plt.xticks(Kpath, x_label)

		lw=1.5
		for i in xrange(len(Ek)):
			ek = array(Ek[i])
			plt.plot(Kpath,ek,"-",color=(0,0.1,0.5),linewidth=lw)

		for kvlp in KVerticalLinePos:
			plt.axvline(x=kvlp, linestyle=":",linewidth=1, color='k')



		plt.show()

if __name__ == "__main__":
	band = Band(sys.argv[1])
	band.plot()
