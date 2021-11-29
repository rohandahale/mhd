# advection.py 

# by Rohan Dahale, November 2021

"""
This is an INPUTS file tailored for the studying. 

NOTE: NX,NY,NZ,DX,DY, and DZ are required constants to inititate the my fluid object; additional problem specific constants are added after creating the my fluid object below.
"""

import numpy as np
import os
cur_path = os.getcwd()

import sys
sys.path.insert(0, '../..')
sys.path.insert(0, '../../..')

from the_fluid_solver_2d import The_Fluid_Solver_2D 
from the_state_saver import save_state

##General Constants
NX= 81 		#number of zones in x
NY= 81	 	#number of zones in y
NZ = "NA"   #not applicable to a 2D problem

XMIN = 0
XMAX = 1.0

YMIN = 0
YMAX = 1.0

DX=(XMAX-XMIN)/(NX-1) #width of a cell in x 
DY=(YMAX-YMIN)/(NY-1) #width of a cell in y 
DZ = "NA"
#End General Constants

my_fluid = The_Fluid_Solver_2D(NX,NY,NZ, DX,DY,DZ)

##Problem Constants
my_fluid.cycle_start = 1                		# start cycle ; use 1 as default 
my_fluid.S 			 = 1E-1             		# stability constant(CFL value)
my_fluid.NI          = 300      			    # number of iterations

my_fluid.RHO         = 1000						# density of the fluid
my_fluid.MU          = 1.0              		# magnetic permeability
my_fluid.NU          = 0                		# viscous diffusion
my_fluid.ETA         = 0                		# resistivity (mag diffusion)

my_fluid.DT          = 5E-3
my_fluid.NT  = 800  	# number of time steps
#End Problem Constants

my_fluid.SAVE_FREQ = int(0.01/my_fluid.DT)  

##Plotting Constants
my_fluid.XMIN = XMIN 
my_fluid.XMAX = XMAX


my_fluid.YMIN = YMIN 
my_fluid.YMAX = YMAX

my_fluid.PLT_TYP = """quiver pcolormesh concentration and mag"""
my_fluid.LABEL = "advection"

my_fluid .MAX_CYC_MAG = 9 # magnitude of max cycles 
##End Plotting Constants

os.chdir(cur_path + "/StateFiles") 
save_state(my_fluid, "velocity.p")

#end advection.py


