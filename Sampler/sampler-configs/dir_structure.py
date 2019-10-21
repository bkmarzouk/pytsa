""" Identifies sampler directory structure """


import os
import sys


# Overide default save location, may want to use scratch on hpc runs
override_save_path = None


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
		os.getcwd(), "..", "sampler-builds", "{}"
	    )
    )
)


# Get default save location for outputs
svepath = (
    os.path.join(
        smppath, "save-files"
    )
)


def get_sampler_paths(build_name):

    if override_save_path is None:
        return pytpath, smppath.format(build_name), svepath.format(build_name)

    return pytpath, smppath.format(build_name), override_save_path
