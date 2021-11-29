import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 18})
import cPickle
import numpy as np
import os
os.chdir("./StateFiles")


the_data = cPickle.load(open("ot_vortex_cycle_0.p","rb"))

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

plt.ion()
plt.contourf(X, Y, B, 30, cmap='RdGy')
plt.colorbar()
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Magnetic Field Intensity at t = 0)')
plt.draw()
plt.pause(1)

Nt = 1000/20
cont = plt.contourf(X, Y, B, 30, cmap='RdGy')


for i in range(20,2000,20):

	t = my_fluid.DT*i

	file = "ot_vortex_cycle_" +str(i) + ".p"
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
	BE = 1.0/2 * B**2 / my_fluid.MU
	U = np.sqrt(u**2 + v**2)

	plt.contourf(X, Y, B, 20, cmap='RdGy')
	plt.quiver(X,Y,u,v)
	plt.xlabel('X')
	plt.ylabel('Y')
	plt.title('Magnetic Field Intensity at t = '+str(t))
	plt.show(block=False)
	plt.pause(0.1)
	plt.gca().cla()



"""
import matplotlib.animation as animation
plt.rcParams['animation.ffmpeg_path'] = '/Volumes/Rohan_1TB/compmhd/'

# animation function
def animate(i):
	global cont
	for c in cont.collections:
		c.remove()  # removes only the contours, leaves the rest intact

	file = "ot_vortex_cycle_" +str(i) + ".p"
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
	U = np.sqrt(u**2 + v**2)

	cont = plt.contourf(X, Y, B, 30, cmap='RdGy')
	return cont
	#plt.quiver(X,Y,u,v)
	#plt.show(block=False)

    #plt.title('t = %i:  %.2f' % (i,z[5,5]))


anim = animation.FuncAnimation(fig, animate, frames=Nt)
"""