#!/usr/bin/env pythonw
from visual import *
from numpy import *
import lift
import sys

cell = lift.Cell()
opt_flag = ""
if		len(sys.argv) == 1:	cell.load('input.lat', 'c')
elif	len(sys.argv) == 2:	cell.load(sys.argv[1], 'o')
elif	len(sys.argv) == 3:
	cell.load(sys.argv[1], 'o')
	opt_flag = sys.argv[2]

print opt_flag


scene=display(title="",x=0,y=0,width=800, height=800, center=(0,0,0), background=(1,1,1))

distant_light(direction=(-1,-1, 0.5), color=(.8,.8,.8))

def get_center(lat):
	max_x = 0
	min_x = 0
	max_y = 0
	min_y = 0
	max_z = 0
	min_z = 0

	for i in xrange(len(lat)):
		atom = cell.superlattice[i]
		UCIndex, AtomIndex, Name_Sub, pos = atom
		x=pos[0,0]
		y=pos[0,1]
		z=pos[0,2]

		if i == 0:
			max_x = x
			min_x = x
			max_y = y
			min_y = y
			max_z = z
			min_z = z
		else:
			if max_x < x : max_x = x
			if max_y < y : max_y = y
			if max_z < z : max_z = z
			if min_x > x : min_x = x
			if min_y > y : min_y = y
			if min_z > z : min_z = z

	#local_light(pos=(-max_x,-max_y, max_z), color=(.8,.8,.8))
	return (min_x,max_x,min_y,max_y,min_z,max_z)

min_x, max_x, min_y, max_y, min_z, max_z = get_center(cell.superlattice)

avg_x = 0.5*(max_x+min_x)
avg_y = 0.5*(max_y+min_y)
avg_z = 0.5*(max_z+min_z)

def plot_cellbox(lat, flag=1):
	radius = 0.02
	A1 = array(cell.a1*cell.supercell_dim[0]).reshape(-1,)
	A2 = array(cell.a2*cell.supercell_dim[1]).reshape(-1,)
	A3 = array(cell.a3*cell.supercell_dim[2]).reshape(-1,)

	x0 = min_x #- avg_x
	y0 = min_y #- avg_y
	z0 = min_z #- avg_z

	d = 0.5
	shift = array([avg_x+d, avg_y+d, avg_z+d])
	p0 = array([x0,y0,z0])

	if flag == 1:
		curve(pos=[p0    - shift, A1   -shift, A1+A2   -shift, A2   -shift, p0   -shift], radius=radius)
		curve(pos=[p0+A3 - shift, A1+A3-shift, A1+A2+A3-shift, A2+A3-shift, p0+A3-shift], radius=radius)

		curve(pos=[p0    - shift, p0+A3-shift]   , radius=radius)
		curve(pos=[A1    - shift, A1+A3-shift]   , radius=radius)
		curve(pos=[A1+A2 - shift, A1+A2+A3-shift], radius=radius)
		curve(pos=[A2    - shift, A2+A3-shift]   , radius=radius)


def cmap(atom_name):
	color_map = {
		"LS": ((0.6,0.6,0.9)	,0.0, 0.5),
		"La": ((0.5,0.5,1.0)	,0.1, 0.5),
		"Sr": ((1.0,0.4,0.0)	,0.1, 0.5),
		"Mn": ((1.0,1.0,0.0)	,0.1, 0.7),
		"Bi": ((0.0,0.6,0.9)	,0.0, 0.5),
		"Fe": ((0.0,1.0,0.0)	,0.1, 1.0),
		"_Fe":((0.0,1.0,0.0)	,0.1, 0.7),
		"O" : ((1.0,1.0,1.0)	,0.0, 0.3),
		"VA": ((0.0,0.0,0.0)	,0.0, 0.1),
		"VC": ((1.0,0.0,0.0)	,0.1, 0.2),
		"C" : ((0.0,0.0,0.0)	,0.1, 0.2),
	}
	return color_map[atom_name]

plot_cellbox(cell.superlattice,1)

radius = 0.12

flag = 1
for atom in cell.superlattice:
	UCIndex, AtomIndex, Name_Sub, pos = atom
	atom_name = Name_Sub.split('-')[0]

	dd = 1
	if len(cell.den_list) > 0:
		dd = sum( cell.den_list[AtomIndex] )

	clr, rad, opa = cmap(atom_name)

	opatune = 1.0

	if pos[0,0] != 4 and pos[0,0] !=4.5 and pos[0,0] !=5 :
	#if pos[0,0] != 4 and pos[0,0] !=4.5 :
		opa = opatune

	location = (pos[0,0]-avg_x, pos[0,1]-avg_y, pos[0,2]-avg_z)
	sphere(pos=location, radius=rad, color=clr, opacity=opa)

	if atom_name == opt_flag:
		label(pos=location, color=color.black, height=10, yoffset=5, text=str(AtomIndex))

	if len(cell.spin_list) > 0:
		if AtomIndex in cell.spin_list.keys():
			sx, sy, sz= cell.spin_list[AtomIndex]

			a = 0.4
			Sx = a*sx.real
			Sy = a*sy.real
			Sz = a*sz.real
		

			opa2 = 0.6
			if pos[0,0] != 4 and pos[0,0] !=4.5 and pos[0,0] !=5 :
			#if pos[0,0] != 4 and pos[0,0] !=4.5 :
				opa2 = opatune

			a = 1.0

			arrow(pos=(pos[0,0]-avg_x, pos[0,1]-avg_y, pos[0,2]-avg_z),
				  axis=(a*Sx,a*Sy,a*Sz), shaftwidth=0.05, color=color.black, opacity=opa2)

			#if atom_name == 'Mn':
			#	arrow(pos=(pos[0,0]-avg_x, pos[0,1]-avg_y, pos[0,2]-avg_z), axis=(a*Sx,0,0),	shaftwidth=0.04, color=(1,0,0), opacity = opa2)
			#	arrow(pos=(pos[0,0]-avg_x, pos[0,1]-avg_y, pos[0,2]-avg_z), axis=(0,a*Sy,0),	shaftwidth=0.04, color=(0,1,0), opacity = opa2)
			#	arrow(pos=(pos[0,0]-avg_x, pos[0,1]-avg_y, pos[0,2]-avg_z), axis=(0,0,a*Sz),	shaftwidth=0.04, color=(0,0,1), opacity = opa2)

			#if atom_name == 'Fe':
			#	arrow(pos=(pos[0,0]-avg_x, pos[0,1]-avg_y, pos[0,2]-avg_z), axis=(a*Sx,0,0),	shaftwidth=0.04, color=(0.7,0,0), opacity = opa2)
			#	arrow(pos=(pos[0,0]-avg_x, pos[0,1]-avg_y, pos[0,2]-avg_z), axis=(0,a*Sy,0),	shaftwidth=0.04, color=(0,0.7,0), opacity = opa2)
			#	arrow(pos=(pos[0,0]-avg_x, pos[0,1]-avg_y, pos[0,2]-avg_z), axis=(0,0,a*Sz),	shaftwidth=0.04, color=(0,0,0.7), opacity = opa2)

