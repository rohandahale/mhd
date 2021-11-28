import matplotlib.pyplot as plt
import cPickle
import numpy as np
import os
os.chdir("./StateFiles")


the_data = cPickle.load(open("ot_vortex_cycle_000000080.p","rb"))

my_fluid = cPickle.load(open("my_fluid.p", "rb"))

X = np.linspace(my_fluid.XMIN,my_fluid.XMAX,my_fluid.NX)
Y = np.linspace(my_fluid.YMIN,my_fluid.YMAX,my_fluid.NY) 
Y,X = np.meshgrid(X,Y)


u 	= the_data[1]
v 	= the_data[2]
p 	= the_data[4]
src = the_data[5]
Bx 	= the_data[8] 
By 	= the_data[9]
A 	= the_data[10]

B = np.sqrt(Bx**2 + By**2)


plt.contourf(X, Y, B, 30, cmap='RdGy')
plt.colorbar()
plt.show(block=False)
plt.pause(0.1)


for i in range(100,900,20):

	file = "ot_vortex_cycle_000000" +str(i) + ".p"
	the_data = cPickle.load(open(file,"rb"))

	my_fluid = cPickle.load(open("my_fluid.p", "rb"))

	X = np.linspace(my_fluid.XMIN,my_fluid.XMAX,my_fluid.NX)
	Y = np.linspace(my_fluid.YMIN,my_fluid.YMAX,my_fluid.NY) 
	Y,X = np.meshgrid(X,Y)


	u 	= the_data[1]
	v 	= the_data[2]
	p 	= the_data[4]
	src = the_data[5]
	Bx 	= the_data[8] 
	By 	= the_data[9]
	A 	= the_data[10]

	B = np.sqrt(Bx**2 + By**2)


	plt.contourf(X, Y, B, 30, cmap='RdGy')
	plt.show(block=False)
	plt.pause(0.1)
