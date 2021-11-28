# the_state_saver.py

# by Rohan Dahale, November 2021

"""
The purpose of this file is to save the state of the fluid processed by the code by pickling the data.

NOTE: using cpickle can be up to 1000x faster than pickle; also note pickle is not secure so do not deserialize from untrusted sources; also note that derserializing won't work if changes to the class are made, so this is not good for long-term storage. 

To load data from a state file, use 

my_fluid = cPickle.load(open("file name.p","rb"))
"""

import cPickle

def save_state(the_data, file_name):
	"""
	Takes as input an array (alone or an array of arrays) and saves it by pickling to file of name 'file_name'.
	"""

	cPickle.dump(the_data, open(file_name, "wb")) 
# end_save_state

# end_the_state_saver.py