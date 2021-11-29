# run_advection.py 

# by Rohan Dahale, November 2021

"""
The purpose of this program is to study advection.
"""

import glob
import numpy as np
import sys
sys.path.insert(0, "../../..")
from the_file_name_getter import get_file_name
from the_state_saver import save_state
from the_fluid import The_Fluid
import os
os.chdir("./StateFiles")
import cPickle
my_fluid = cPickle.load(open("velocity.p", "rb"))

##Setup
ext = ".p" # filename extension for pickling

X = np.linspace(my_fluid.XMIN,my_fluid.XMAX,my_fluid.NX)
lx=len(X)
Y = np.linspace(my_fluid.YMIN,my_fluid.YMAX,my_fluid.NY) 
Y,X = np.meshgrid(X,Y)

XX,YY= np.mgrid[1:my_fluid.NX-1, 1:my_fluid.NY-1] #imitates 2d for loop

if my_fluid.cycle_start == 1:
	u 	= np.zeros((my_fluid.NX,my_fluid.NY)) # x-component of velocity
	v 	= np.zeros((my_fluid.NX,my_fluid.NY)) # y-component of velocity
	p 	= np.zeros((my_fluid.NX,my_fluid.NY)) # pressure
	src = np.zeros((my_fluid.NX,my_fluid.NY)) # src term for poisson eqn
	Bx 	= np.zeros((my_fluid.NX,my_fluid.NY)) # x-comp of magnetic field
	By 	= np.zeros((my_fluid.NX,my_fluid.NY)) # y-comp of magnetic field
	A 	= np.ones((my_fluid.NX,my_fluid.NY)) # magnetic vector potential

	##Initialize single vortex for velocity; double vortex for mag field

	
	u = np.cos(0*np.pi*Y)
	v = np.sin(0*np.pi*X)

	Bx 	= 2*np.sin(0*np.pi*X)
	By 	= 2*np.cos(0*np.pi*Y)

	B0 = 1.0
	A  = (X + Y)/np.sqrt(2)

	#Update ghost zones so boundaries are periodic
	u[ 0, : ] 	= u[ -2, : ]
	u[ -1, : ] 	= u[ 1, : ] 
	u[ :, 0 ]	= u[ :, -2 ]
	u[ :, -1 ] 	= u[ :, 1 ]

	v[ 0, : ] 	= v[ -2, : ]
	v[ -1, : ] 	= v[ 1, : ]
	v[ :, 0 ] 	= v[ :, -2 ]
	v[ :, -1 ] 	= v[ :, 1 ]

	A[ 0, : ] 	= A[ -2, : ]
	A[ -1, : ] 	= A[ 1, : ]
	A[ :, 0 ] 	= A[ :, -2 ]
	A[ :, -1 ]  = A[ :, 1 ]

	Bx[1: -1 ,1: -1] = (A[1: -1 ,2:] - A[1: -1 ,: -2])/my_fluid.DY 
	By[1: -1 ,1: -1] = -(A[2: ,1: -1] - A[: -2 ,1: -1])/my_fluid.DX

	#Update ghost zones so boundaries are periodic
	Bx[ 0, : ] 	= Bx[ -2, : ]
	Bx[ -1, : ] = Bx[ 1, : ]
	Bx[ :, 0 ]  = Bx[ :, -2 ]
	Bx[ :, -1 ] = Bx[ :, 1 ]

	By[ 0, : ]  = By[ -2, : ]
	By[ -1, : ] = By[ 1, : ] 
	By[ :, 0 ]  = By[ :, -2 ]
	By[ :, -1 ] = By[ :, 1 ]

	##Save the state of the initial condition of the fluid
	cycles = 0
	file_name = get_file_name(cycles, my_fluid.MAX_CYC_MAG, my_fluid.LABEL, ext)

	save_state([cycles, u,v, "NA", p,src, "NA","NA", Bx,By,A], file_name)

else:
	data_file = glob.glob("*" + str(my_fluid.cycle_start) + "*")

	if data_file == []:
		print("Error! my fluid.cycle_start data file not found.")
		sys.exit()
	#end if
	data_file = data_file[0]
	the_data  = cPickle.load(open(data_file))

	u 	= the_data[1]
	v 	= the_data[2]
	p 	= the_data[4]
	src = the_data[5]
	Bx 	= the_data[8] 
	By 	= the_data[9]
	A 	= the_data[10]
	the_data = 0 #to save memory
#end if


#Solve

for cycles in xrange(my_fluid.cycle_start, my_fluid.NT+1): 
	u_old = u.copy()

	u[1:-1,1:-1] = (my_fluid.nonlinear_advect_implicit_periodic_2d(u,u,v, XX,YY) - my_fluid.nonlinear_advect_implicit_periodic_2d(Bx,Bx,By, XX,YY) + Bx[1: -1 ,1: -1] - my_fluid.apply_pressure_2dX(Bx**2+By**2, 0.5))

	v[1:-1,1:-1] = (my_fluid.nonlinear_advect_implicit_periodic_2d(v,u_old,v, XX,YY) - my_fluid.nonlinear_advect_implicit_periodic_2d(By,Bx,By, XX,YY) + By[1: -1 ,1: -1] - my_fluid.apply_pressure_2dY(Bx**2+By**2, 0.5))

	##Update ghost zones
	u[ 0, : ]  = u[ -2, : ]
	u[ -1, : ] = u[ 1, : ]
	u[ :, 0 ]  = u[ :, -2 ]
	u[ :, -1 ] = u[ :, 1 ]

	v[ 0, : ]  = v[ -2, : ]
	v[ -1, : ] = v[ 1, : ]
	v[ :, 0 ]  = v[ :, -2 ]
	v[ :, -1 ] = v[ :, 1 ]
	src[1:-1, 1:-1] = my_fluid.calc_source_2d(u,v)
	p = my_fluid.transform_pressure_poisson_2d(p, src)

	p[ 0,:]	  = p[ -2, :]
	p[ -1, :] = p[ 1,: ]
	p[ :, 0 ] = p[ :, -2 ]
	p[ :,-1 ] = p[ :, 1 ]
 
	u[1:-1,1:-1] -= my_fluid.apply_pressure_2dX(p, 1.0/my_fluid.RHO)
	v[1:-1,1:-1] -= my_fluid.apply_pressure_2dY(p, 1.0/my_fluid.RHO)

	##Update ghost zones
	u[ 0, : ]  = u[ -2, :]
	u[ -1, : ] = u[ 1, : ]
	u[ :, 0 ]  = u[ :, -2]
	u[ :, -1]  = u[ :, 1 ]

	v[ 0, : ]  = v[ -2, :]
	v[ -1, : ] = v[ 1, : ]
	v[ :, 0 ]  = v[ :, -2]
	v[ :, -1]  = v[ :, 1 ]

	A[1: -1 ,1: -1] = my_fluid.nonlinear_advect_implicit_periodic_2d(A,u,v, XX,YY)

	A[ 0, : ]  = A[ -2, : ]
	A[ -1, : ] = A[ 1, : ] 
	A[ :, 0 ]  = A[ :, -2 ]
	A[ :, -1 ] = A[ :, 1 ]


	Bx[1:-1,1:-1] = ((A[1:-1, 2: ] - A[1:-1, :-2])/my_fluid.DY)
	By[1:-1,1:-1] = -((A[2: , 1:-1] - A[ :-2, 1:-1])/my_fluid.DX)

	Bx[ 0, : ]  = Bx[ -2, : ]
	Bx[ -1, : ] = Bx[ 1, : ]
	Bx[ : 0 ]   = Bx[ :, -2 ]
	Bx[ :, -1 ] = Bx[ :, 1 ]

	By[ 0, : ]  = By[ -2, : ]
	By[ -1, : ] = By[ 1, : ]
	By[ : 0 ]   = By[ :, -2 ]
	By[ :, -1 ] = By[ :, 1 ]

	if (cycles%my_fluid.SAVE_FREQ == 0):
		file_name = get_file_name(cycles, my_fluid.MAX_CYC_MAG, my_fluid.LABEL, ext) 
		save_state([cycles, u,v, "NA", p,src, "NA","NA", Bx,By, A],file_name)
	#end if
#end for

file_name = get_file_name(cycles, my_fluid.MAX_CYC_MAG, my_fluid.LABEL, ext)
save_state([cycles, u,v, "NA", p,src, "NA","NA", Bx,By,A], file_name)

# end implicit_orszag_tang_athena_2d.py



