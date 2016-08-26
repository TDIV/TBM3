![starscream](https://raw.githubusercontent.com/TDIV/TBM3/master/icon/icon.png)

TBM3 (Tight Binding Model for Materials at Mesoscale) is a C++ based numerical tool package and framework that designed for construct ‘any’ kind of lattice structure with multi-orbital and spin degree of freedom to solve the basic nano-scale electronic structure. It is to designed to construct a highly flexible and reusable framework for anyone to easily solve the multi-scale quantum mechnical model and extract useful physical quantities.

## Features

- A general lattice input file with orbital information specified.
- Providing non-spin polarized and spin polarized calculation with normal (non-superconductive), Nambu and extended Nambu space.
- Meanfield self-consistence and LLG spin dynamic calculations.
- Simple and minimal input format for the setting of a quantum mechanical model.
- Simple file formate for the storage of order parameters.

## Simple Example

First thing is to construct the lattice file. Here is an example for a single site Fe atom with two orbitals (dxz & dyz) for the iron based compound (BaFe2As2). Here we name the lattice input file as `BaFe2As2.lat`.

```
#BasisVector
 1               0               0              
 0               1               0              
 0               0               1              

#OrbitalProfile
Fe        dxz       dyz       

#Atoms
1     0               0               0      
```

Next, we can design and construct the tight-binding hopping term for the lattice. Here we choose the parameter sets from this paper [Phys. Rev. B 88, 184509 (2013)](http://journals.aps.org/prb/abstract/10.1103/PhysRevB.88.184509). Here we the model setting input file follows the lattice filename as `BaFe2As2.lat.tbm`.

```
#Parameters
isCalculateMu	= 1
isCalculateVar 	= 0
isCalculateLDOS	= 0
isCalculateBand	= 1

spin		= "on"
space		= "normal"
bondRadius	= 3

Mu			= 0
Temperature	= 0.0001

t1			= 0.09
t2			= 0.08
t3			= 1.35
t4			= -0.12
t5			= -1
t6			= 0.25

#CoreCharge
Fe	> 2

#Hamiltonian
hoppingHc > Fe:Fe:+1+0+0 1:1 >  t1
hoppingHc > Fe:Fe:+1+0+0 2:2 >  t1
hoppingHc > Fe:Fe:+0+1+0 1:1 >  t1
hoppingHc > Fe:Fe:+0+1+0 2:2 >  t1

hoppingHc > Fe:Fe:+1+1+0 1:1 >  t2
hoppingHc > Fe:Fe:+1-1+0 2:2 >  t2

hoppingHc > Fe:Fe:+1+1+0 2:2 >  t3
hoppingHc > Fe:Fe:+1-1+0 1:1 >  t3

hoppingHc > Fe:Fe:+1+1+0 1:2 >  t4
hoppingHc > Fe:Fe:+1-1+0 1:2 >  t4
hoppingHc > Fe:Fe:-1-1+0 1:2 >  t4
hoppingHc > Fe:Fe:-1+1+0 1:2 >  t4

hoppingHc > Fe:Fe:+1+0+0 1:2 >  t5
hoppingHc > Fe:Fe:+0+1+0 1:2 >  t5
hoppingHc > Fe:Fe:-1+0+0 1:2 >  t5
hoppingHc > Fe:Fe:+0-1+0 1:2 >  t5

hoppingHc > Fe:Fe:+2+0+0 1:1 >  t6
hoppingHc > Fe:Fe:+0+2+0 1:1 >  t6
hoppingHc > Fe:Fe:+2+0+0 2:2 >  t6
hoppingHc > Fe:Fe:+0+2+0 2:2 >  t6

#KPointPath
G      0     0     0    
X      1     0     0    
M      1     1     0    
G      0     0     0    
```

In above, the keyword `hopingHc` describe a Hermitian conjugate hopping term from one Fe site to another (`Fe:Fe:+1+0+0`) for the specific (inter-)orbitals (`1:1`, `2:2`, `1:2`) with a given parameter (`t1`, `t2`, `t3`, `t4`, `t5`). The bond direction is described in the notation `+1+0+0`.

Afterwards, you can create another block to calculate the band structure along the given `#KPointPath`.

Now, we are good to execute the program under the terminal:

```
tbm-run BaFe2As2.lat
```

It will automatically calculate the chemical potential for the given electron charge (2 electrons per Fe site), and output the band structure 'BaFe2As2.lat.ban' under the same folder.

## Requirements
- A NVIDIA GPU.
- CUDA library.
- MAGMA library.
- BOOST library.

We are working on to provide more possible options for different hardward selections.

## Installation

Some installation instructions.

## License

TBM3 is licensed under the GNU General Public License.

## Contact

- Dr. Yuan-Yen Tai: [Personal web](http://dr-tai.net), [Linked In](https://www.linkedin.com/in/yuan-yen-tai-5652ab112)

