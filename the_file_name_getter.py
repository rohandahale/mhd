# the_file_name_getter.py

# by Rohan Dahale, November 2021

"""
Prepends zeros to the cycle number to be used in the file saved name so that files show up in the correct order on drive.
"""

def get_file_name(cycles, max_cyc_mag, label, ext):
	"""
	Requires as input the current cycle (as float) and the maximum number of cycles (as float) that will be run; also requires as input the my fluid 'LABEL' attribute (as string) as well as the file name extension (as string).

	Ouputs the cycle label as a string.
	"""

	cyc_lbl = ""
	for i in xrange(1, max_cyc_mag + 1):
		if cycles < 10**i:
			#cyc_lbl = "0" *(max_cyc_mag - i) + str(cycles)
			cyc_lbl =  str(cycles)
			file_name = label + "_cycle_" + cyc_lbl + ext
			return file_name 
		# end if
	# end for
# end get_file_name

# end the_file_name_getter.py