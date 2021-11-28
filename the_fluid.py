#the_fluid.py

# By Rohan Dahale, Novemeber 2021


'''
The Fluid contains all the instance wide constants (attributes) that will be used by this code to solve fluid mechanics problems; these attributes and The_Fluid instance are instantiated by the input file of a problem.
'''


class The_Fluid(object):
	def __init__(self, NX,NY,NZ, DX,DY, DZ):
		self.NX = NX # number of zones in x
		self.NY = NY
		self.NZ = NZ
	
		self.DX = DX # width of a cell(zone) in x
		self.DY = DY
		self.DZ = DZ
	# end __init__
#end class The_Fluid

# end the_fluid.py
