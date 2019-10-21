""" Constructs sampler directory and programme """


# Import system modules
import sys
import os
import shutil
import importlib

# Get python path
pytpath = (
    os.path.abspath(
        os.path.join(
            os.getcwd(), "../..", "PyTransport", "PyTransport"
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

sys.path.append(pytpath)
sys.path.append(smppath)


# Setup PyTransport installation paths
import PyTransSetup as PySet
PySet.pathSet()


# Load user arguments
user_args = sys.argv
assert len(user_args) == 2, "command: 'python build_sampler.py <config>.py'"


# Get configuration module of interest
config_module = user_args[1]


# Import as module based on cmd line arg
config = importlib.import_module(config_module[:-3])


# Check that module is installed
sampler = config.sampler
try:
    PyT = importlib.import_module(sampler["PyTransportModule"])

except ImportError:
    raise ImportError, "Failed to import PyTransport module: {}, check installation.".format(
        sampler["PyTransportModule"]
    )


# Get required modules for sample generator
modules = config.required_modules


# Load file ICs and model parameters
field_positions  = config.field_positions
field_velocities = config.field_velocities
parameter_values = config.parameter_values


# Get installation number for cross validation
nF = PyT.nF()
nV = PyT.nF()
nP = PyT.nP()


# Build set of fiducial values for fields, velocities and parameters
fiducial_fvals = [None for i in range(nF)]; cf = 0
fiducial_vvals = [None for i in range(nV)]; cv = 0
fiducial_pvals = [None for i in range(nP)]; cp = 0


# Read each field value: if "ALL" is found, replace fiducial values with method prescribed in confiuration
for p in field_positions:
    if p["FieldNumber"] == "ALL":
        fiducial_fvals = [p["Command"] for i in range(nF)]


# Now reiterate over items: If anything but "ALL" is found, replace the fiducial method for the field
for p in field_positions:
    f_num = p["FieldNumber"]
    if f_num != "ALL":
        fiducial_fvals[f_num] = p["Command"]
assert None not in fiducial_fvals


# Repeat method for velocities
for v in field_velocities:
    if v["FieldNumber"] == "ALL":
        fiducial_vvals = [v["Command"] for i in range(nV)]
for v in field_velocities:
    f_num = v["FieldNumber"]
    if f_num != "ALL":
        fiducial_vvals[f_num] = v["Command"]
assert None not in fiducial_vvals


# Repeat method for paramaters
for p in parameter_values:
    if p["ParameterNumber"] == "ALL":
        fiducial_pvals = [p["Command"] for i in range(nP)]
for p in parameter_values:
    p_num = p["ParameterNumber"]
    if p_num != "ALL":
        fiducial_pvals[p_num] = p["Command"]
assert None not in fiducial_pvals


# We construct strings for generator files: We ned the correct formatting to be returned by the generator
f_string = "["; fc = 1
v_string = "["; vc = 1
p_string = "["; pc = 1


for f in fiducial_fvals:
    if type(f) is not str:
        f_string += str(f)
    else: f_string += f
    if fc != nF:
        fc+=1
        f_string += ","
    else: f_string += "]"


for v in fiducial_vvals:
    if type(v) is not str:
        v_string += str(v)
    elif v == "SlowRoll": v_string+="'{}'".format(v)
    else: v_string += v
    if vc != nV:
        vc+=1
        v_string += ","
    else: v_string += "]"


for p in fiducial_pvals:
    if type(p) is not str:
        p_string += str(p)
    else: p_string += p
    if pc != nP:
        pc+=1
        p_string += ","
    else: p_string += "]"


# Setup import statements for writing process
import_lines   = ["import {}\n".format(m) for m in modules]
function_lines = [
    "def gen_sample(n):\n",
    "\tfvalues = {}\n".format(f_string),
    "\tvvalues = {}\n".format(v_string),
    "\tpvalues = {}\n".format(p_string),
    "\treturn n, fvalues, vvalues, pvalues\n"
]


# Build save directory
saveloc = config.save_location
sampler = config.sampler_name

if saveloc == "default":
    assert os.path.exists(smppath), "sampler-build directory not found: {}".format(smppath)
    savepath = os.path.join(smppath, sampler)
else:
    assert os.path.exists(saveloc), "directory for save files not found: {}".format(saveloc)
    savepath = os.path.join(saveloc, sampler)

    
assert not os.path.exists(savepath), "Directory already exists! {}".format(saveloc)
os.makedirs(savepath)


# Write sample generator file
gfile = open(os.path.join(savepath, "generator.py"), "w")
with gfile:
    print "-- Writing sample generator"
    for line in import_lines + function_lines:
        gfile.writelines(line)
        
# Write configuration file, with new directory information
cfg_new_file = open(os.path.join(savepath, "config.py"), "w")
cfg_old_file = open(config_module, "r")

cfg_lines = []

with cfg_old_file as c:
    for line in c.readlines():
        if "sys01" in line:
            line_ = "\t'pytpath': '{}',\n".format(pytpath)
            cfg_lines.append(line_)
            print line_
        elif "sys02" in line:
            line_ = "\t'saveloc': '{}',\n".format(savepath)
            cfg_lines.append(line_)
            print line_
        elif "sys03" in line:
            line_ = "\t'smppath': '{}'".format(os.path.abspath(os.path.join(pytpath, "../../Sampler")))
            cfg_lines.append(line_)
            print line_
        elif "save_location" in line:
            pass
        elif "sampler_name" in line:
            pass
        else:
            cfg_lines.append(line)
            
with cfg_new_file as c:
    for line in cfg_lines:
        c.writelines(line)


# Locate pre-written sampler files: e.g. mapreduce, etc
sampler_file_dir = os.path.abspath(os.path.join(
    os.getcwd(), "..", "sampler-methods"
))

sampler_files = [
    "run_sampler.py"
]


# Copy appropriate files
for item in sampler_files:
    print "-- Copying {}".format(item)
    shutil.copyfile(
        os.path.join(sampler_file_dir, item), os.path.join(savepath, item)
    )
# shutil.copyfile(os.path.abspath(os.path.join(os.getcwd(), config_module)),
#                 os.path.join(saveloc, "config.py"))



# Construct ensemble of
build_dirs = ["samples", "outputs"]
if config.computations['2pf'] is True: # 2-point function data
    build_dirs.append("2pf")

if config.computations['3pf'] is True: # 3-point function data
    build_dirs.append("3pf")

if config.computations['Mij'] is True: # Mass-Matrix data, (at Horizon crossing)
    build_dirs.append("Mij")


for item in build_dirs:

    d = os.path.join(savepath, item)

    assert not os.path.exists(d)
    
    os.makedirs(d)
