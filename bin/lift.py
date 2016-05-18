#!/usr/bin/env python

import numpy as np
from math import *
import coord_search as cs

# -------------------------------------------------------------
# ---String operations-----------------------------------------
# -------------------------------------------------------------
def tos(var, length=10):
	var = str(var)
	if var[0] != '-' and var[0] != '+':
		var = ' ' + var
	while len(var) < length - 1:
		var += ' '
	return var + ' '

# Formate the line with specific length
def line_formate	(line, length_list):
	line=line.split()
	ret_line=""
	for ii in range(len(line)):
		if type(length_list) is list:
			ret_line += tos(line[ii], length_list[ii])
		else:
			ret_line += tos(line[ii], length_list)
	return ret_line
def parse_bond_str	(bonding_basis, word):
	bond_vec = np.mat([0, 0, 0])
	#word = word.replace(".", "+0")
	tmpWord = word[1:]
	tmpWord = tmpWord.replace("+",",")
	tmpWord = tmpWord.replace("-",",")
	Nxyz = tmpWord.split(',')
	Nx = int(Nxyz[0])
	Ny = int(Nxyz[1])
	Nz = int(Nxyz[2])
	signArray = []
	for s in word:
		if s == '+':
			signArray.append(+1)
		elif s == '-':
			signArray.append(-1)

	#print bonding_basis[0][1]
	bond_vec += signArray[0] * Nx * bonding_basis[0][1]
	bond_vec += signArray[1] * Ny * bonding_basis[1][1]
	bond_vec += signArray[2] * Nz * bonding_basis[2][1]
	#print bond_vec, word

	return bond_vec
def cpp_form		(val):
	if isinstance(val, int):
		return str(val)
	elif isinstance(val, float):
		return str(val)
	elif isinstance(val, complex):
		return "(" + str(val.real) + "," + str(val.imag) + ")"
	else:
		return str(val)

def parse_order( word ):
	sub  = word.split(':')
	sval = sub[1][1:-1]
	sp_ang = sval.split(',')
	cval= float(sp_ang[0])+1j*float(sp_ang[1])
	return (sub[0], cval)
# -------------------------------------------------------------
# The class for parsing the input parameters-------------------
# -------------------------------------------------------------
class ParseParameter:

	def __init__(self):
		self.parameterList = []
		self.parameter = {}

	def parse(self, line):
		line = line.replace(" ", "")
		line = line.replace("\t", "")
		line = line.replace("\n", "")

		ll = line.split("=")

		if len(ll) == 2:

			if len(self.parameterList) == 0:
				tmp_value = eval(ll[1])
				self.parameter[ll[0]] = tmp_value
				self.parameterList.append((ll[0], str(tmp_value)))
			else:
				for name, val in self.parameterList:
					ll[1] = ll[1].replace(name, str(val))

				tmp_value = eval(ll[1])

				self.parameter[ll[0]] = tmp_value
				self.parameterList.append((ll[0], str(tmp_value)))

			return tos(ll[0], 16) + "=" + \
				tos(cpp_form(self.parameter[ll[0]]), 16)

		if len(ll) == 1 > 0:

			if len(ll[0]) > 0:
				for name, val in self.parameterList:
					ll[0] = ll[0].replace(name, str(val))

				tmp_value = eval(ll[0])
				return str(tmp_value)

		return ""

	def eval(self, line):
		return eval(self.parse(line))

# -------------------------------------------------------------
# ---Construct the lattice class-------------------------------
# -------------------------------------------------------------
class Cell:

	def __init__(self):  # Initially read in the file structure
		self.filename = ""
		self.save_name = ""
		self.parameter = []
		self.a1 = []
		self.a2 = []
		self.a3 = []
		self.sub_atom = []
		self.supercell_dim = []
		self.bonding_basis = []
		self.bondings = []
		self.superlattice = []
		self.bond = []
		self.Neighbors = []
		self.formated_lines = []
		self.parser = ParseParameter()

		self.den_list = {}
		self.spin_list = {}

		self.maxNeighbor = 1

	def load(self, filename, file_flag='c'):  # Initially read in the file structure

		self.filename = filename
		self.save_name = filename+'.lif'
		formated_lines = []
		flag = -1
		##################################
		# Start analyze the file
		##################################
		f = open(self.filename)
		for l in f.readlines():

			if l[:2] == '#0': flag = 0; formated_lines.append(l.replace('\n', '')); continue;
			if l[:2] == '#1': flag = 1; formated_lines.append(l.replace('\n', '')); continue;
			if l[:2] == '#2': flag = 2; formated_lines.append(l.replace('\n', '')); continue;
			if l[:2] == '#3': flag = 3; formated_lines.append(l.replace('\n', '')); continue;
			if l[:2] == '#4': flag = 4; formated_lines.append(l.replace('\n', '')); continue;
			if l[:2] == '#5': flag = 5; formated_lines.append(l.replace('\n', '')); continue;
			if l[:2] == '#6': flag = 6; continue;
			if l[:2] == '#7': flag = 7; continue;
			if l[:2] == '#8': flag = 8; continue;

			if flag == 0:
				formated_lines.append(l.replace('\n', ''))
				self.parameter.append(self.parser.parse(l))

			# Read the basis vector
			if flag == 1:
				ss = []
				if l[0] != '#':
					formated_lines.append(line_formate(l, [10, 20, 20, 20]))
					ss = l.split()
				if len(ss) > 0:
					if ss[0] == 'a1': self.a1 = [self.parser.eval(ss[1]), self.parser.eval(ss[2]), self.parser.eval(ss[3])]
					if ss[0] == 'a2': self.a2 = [self.parser.eval(ss[1]), self.parser.eval(ss[2]), self.parser.eval(ss[3])]
					if ss[0] == 'a3': self.a3 = [self.parser.eval(ss[1]), self.parser.eval(ss[2]), self.parser.eval(ss[3])]

			# Read the sub atom structure
			if flag == 2:
				ss = []
				if l[0] != '#':
					formated_lines.append(line_formate(l, [5, 5, 20, 20, 20]))
					ss = l.split()
				if len(ss) > 0:
					self.sub_atom.append((ss[0], ss[1], np.mat(
						[self.parser.eval(ss[2]), self.parser.eval(ss[3]), self.parser.eval(ss[4])])))

			# Read the Super cell dim
			if flag == 3:
				ss = []
				if l[0] != '#':
					formated_lines.append(line_formate(l, [5, 5, 5]))
					ss = l.split()
				if len(ss) > 0:
					self.supercell_dim = [int(ss[0]), int(ss[1]), int(ss[2])]

			# Read the bonding basis
			if flag == 4:
				ss = []
				if l[0] != '#':
					formated_lines.append(line_formate(l, [10, 20, 20, 20]))
					ss = l.split()
				if len(ss) > 0:
					self.bonding_basis.append(
						(ss[0], np.mat([self.parser.eval(ss[1]), self.parser.eval(ss[2]), self.parser.eval(ss[3])])))

			# Read the supercell lattice
			if flag == 7:
				ss= []
				if l[0] != '#':
					#formated_lines.append(line_formate(l, [10, 10, 10, 20, 20, 20]))
					ss = l.split()
				if len(ss) > 0:
					UCIndex, AtomIndex, Name_Sub, px, py, pz = \
						int(ss[0]), int(ss[1]), ss[2], float(ss[3]), float(ss[4]), float(ss[5])
					pos = np.mat([px,py,pz])
					self.superlattice.append(
						(UCIndex, AtomIndex, Name_Sub, pos)
					)

			if flag == 8 and l[0] != '#':
				self.Neighbors.append(l.replace('\n', ''))

		self.formated_lines = formated_lines

		self.a1 = np.mat(self.a1)
		self.a2 = np.mat(self.a2)
		self.a3 = np.mat(self.a3)

		self.generateBondings()

		f.close()

		#self.load_den(self.filename+".den")
		self.load_spin(self.filename+".ord")
		##################################
		# Close the file
		##################################

		if len(self.superlattice) == 0 and len(self.Neighbors) ==0 :
			self.CreateLattice()
			self.CreateNeighbor()
		elif file_flag == 'c':
			self.CreateLattice()
			self.CreateNeighbor()
		elif file_flag == 'cn':
			self.CreateNeighbor()

		self.save()
	def generateBondings(self):
		Nx = self.maxNeighbor
		Ny = self.maxNeighbor
		Nz = self.maxNeighbor
		if np.sum(abs(self.bonding_basis[0][1])) == 0 :
			Nx = 0
		if np.sum(abs(self.bonding_basis[1][1])) == 0 :
			Ny = 0
		if np.sum(abs(self.bonding_basis[2][1])) == 0 :
			Nz = 0

		def NumStr(val):
			if val == 0:
				return "+0"
			elif val > 0:
				return "+"+str(abs(val))
			elif val < 0:
				return "-"+str(abs(val))

		for ix in xrange(-Nx,Nx+1):
			for iy in xrange(-Ny,Ny+1):
				for iz in xrange(-Nz,Ny+1):
					bondStr = NumStr(ix)+NumStr(iy)+NumStr(iz)
					self.bondings.append(bondStr)

	def load_den(self, filename):
		try:
			f = open(filename)

			for l in f.readlines():
				if l[0] != '#':
					ss = l.split()
					orbital_den= []
					if len(ss) > 8:
						orbital_str = ss[8:]
						for ostr in orbital_str:
							val = ostr.split(':')
							orbital_den.append(float(val[1]))
					self.den_list[int(ss[0])] = orbital_den

			f.close()
		except IOError:
			print "Warning: Could not open file: "+filename

	def load_spin(self, filename):
		try:
			f = open(filename)

			for l in f.readlines():
				if l[0] != '#':
					ss = l.split()
					cang = 0
					if len(ss) == 12:
						sx, vx = parse_order(ss[9])
						sy, vy = parse_order(ss[10])
						sz, vz = parse_order(ss[11])
						self.spin_list[int(ss[1])] = (vx,vy,vz)

			f.close()
		except IOError:
			print "Warning: Could not open file: "+filename

	def CreateLattice(self):  # Create a super lattice structure
		UCIndex = 0
		AtomIndex = 0

		dim = self.supercell_dim
		n1, n2, n3 = dim[0], dim[1], dim[2]

		self.superlattice = []
		for i3 in xrange(n3):
			for i2 in xrange(n2):
				for i1 in xrange(n1):
					A1i = self.a1 * i1
					A2i = self.a2 * i2
					A3i = self.a3 * i3

					A123 = A1i + A2i + A3i

					SubAtom_index = 0
					for atom in self.sub_atom:
						Name, orb, p = atom

						pos = p + A123
						self.superlattice.append(
							(UCIndex, AtomIndex, Name + '-' + str(SubAtom_index), pos))
						AtomIndex += 1
						SubAtom_index += 1
					UCIndex += 1

	def CreateNeighbor(self):  # Create the table of index relations to the neighbors
		self.bond = []
		for bd in self.bondings:
			self.bond.append((bd, parse_bond_str(self.bonding_basis, bd)))

		N1 = 2
		N2 = 2
		N3 = 2

		if sqrt(self.a1*self.a1.transpose()) < 0.001: N1=0
		if sqrt(self.a2*self.a2.transpose()) < 0.001: N2=0
		if sqrt(self.a3*self.a3.transpose()) < 0.001: N3=0

		# -------------------------------------------------------------------------
		# Scan through superlattice to create a map[pos]=(UCIndex, Index, Label)
		# Periodicity of the super lattice is considered here
		# -------------------------------------------------------------------------
		cs_superlattice = cs.coord_search()
		for atom in self.superlattice:
			UCIndex, Index, Label, pos = atom

			A1 = self.a1 * self.supercell_dim[0]
			A2 = self.a2 * self.supercell_dim[1]
			A3 = self.a3 * self.supercell_dim[2]
			for n1 in range(-N1, N1 + 1):
				for n2 in range(-N2, N2 + 1):
					for n3 in range(-N3, N3 + 1):
						nn = A1 * n1 + A2 * n2 + A3 * n3
						cs_superlattice.add_point(pos+nn, Index)
		# #########################################################################

		# -------------------------------------------------------------------------
		# Scan through superlattice and construct bonding index
		# -------------------------------------------------------------------------
		self.Neighbors = []
		for atom in self.superlattice:
			UCIndex, Index, Label, pos = atom

			line = tos(Index)

			for b_str, b_vec in self.bond:
				near_pos = pos + b_vec

				index, found = cs_superlattice.search(near_pos)

				if index > -1:
					line += tos(b_str + ":" + str(found.data))

			self.Neighbors.append(line)
		# #########################################################################

	def save(self):

		f_lat = open(self.filename, 'w')
		for l in self.formated_lines:
			f_lat.write(l+'\n')

		f = open(self.save_name, 'w')
		a1 = self.a1
		a2 = self.a2
		a3 = self.a3
		sa1 = tos('a1') + tos(a1[0, 0], 20) + tos(a1[0, 1], 20) + tos(a1[0, 2], 20)
		sa2 = tos('a2') + tos(a2[0, 0], 20) + tos(a2[0, 1], 20) + tos(a2[0, 2], 20)
		sa3 = tos('a3') + tos(a3[0, 0], 20) + tos(a3[0, 1], 20) + tos(a3[0, 2], 20)


		f.write('#0 Parameters\n')
		for sa in self.parameter:
			f.write(sa+'\n')

		f.write('#1 Basis Vector\n')
		f.write(sa1 + '\n')
		f.write(sa2 + '\n')
		f.write(sa3 + '\n')

		f.write('\n')
		f.write('#2 Sub Atoms\n')
		for sa in self.sub_atom:
			name, orb, pos = sa
			sspos = tos(name, 5) + tos(orb, 5) + \
				tos(pos[0, 0], 20) + tos(pos[0, 1], 20) + tos(pos[0, 2], 20)
			f.write(sspos + '\n')

		dim = self.supercell_dim
		f.write('\n')
		f.write('#3 Cell Dim\n')
		f.write(tos(dim[0], 5) + tos(dim[1], 5) + tos(dim[2], 5) + '\n')

		f.write('\n')
		f.write('#4 Bonding basis\n')
		for sa in self.bonding_basis:
			name, pos = sa
			sspos = tos(name) + tos(pos[0, 0], 20) + tos(pos[0, 1], 20) + tos(pos[0, 2], 20)
			f.write(sspos + '\n')

		f.write('\n')
		f.write('#7 Cell\n')
		f_lat.write('#7 Cell\n')
		for sl in self.superlattice:
			AtomIndex, UCIndex, Name, pos = sl
			line = tos(AtomIndex) + tos(UCIndex) + tos(Name) + \
					tos(pos[0, 0], 20) + tos(pos[0, 1], 20) + tos(pos[0, 2], 20) + '\n'
			f.write(line)
			f_lat.write(line)

		f.write('\n')
		f.write('#8 Neighbors\n')
		f_lat.write('\n')
		f_lat.write('#8 Neighbors\n')
		for line in self.Neighbors:
			f.write(line_formate(line, 10) + '\n')
			f_lat.write(line_formate(line, 10) + '\n')

		f.close()
		f_lat.close()

	def getAtomName(self, index):
		satom = self.sub_atom[index]
		return satom[0]

if __name__ == '__main__':

	import sys

	cell = Cell()
	if len(sys.argv) == 1:
		cell.load('input.lat', 'c')

	elif len(sys.argv) == 2:
		cell.load(sys.argv[1], 'o')

	elif len(sys.argv) == 3:
		cell.load(sys.argv[1], sys.argv[2])












