""" Identifies sampler directory structure """

def get_pyt_paths():

	import os
	import sys

	# Get python path
	pytpath = (
	    os.path.abspath(
		    os.path.join(
			os.getcwd(), "..", "PyTransport", "PyTransport"
		    )
	    )
	)


	# Get default location for samplers
	smppath = (
	    os.path.abspath(
		    os.path.join(
			os.getcwd(), "..", "sampler-builds"
		    )
	    )
	)

	return pytpath, smppath
