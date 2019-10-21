"""
    Builds N-flation type models with canonical field space and quadratic/chaotic potential terms

    execute "python NQuadGenerator.py <int>"

    and a <int>:QuadSetup.py file will be generated.

    Executing "python <int>QuadSetup.py" will then construct the PyTransport module.

"""


# We insert the preamble for the potential here
import sys


# Load integer user argument
user = sys.argv
assert len(user) == 2, "Enter the number of field for the N-Quadratic field potential."


# Try to convert to integer
try:
    N = int(user[1])
except:
    print "user arg must be integere, i.e. 'python NQuadTemp.py <number>'"


# Define file name for compilation
fname = str(N) + "QuadSetup.py"


# Now we define al the lines we will write to file
preamble = ["import os\n",
            "import sympy as sym\n",
            "import sys\n",
            "tree = os.path.abspath(os.path.join(__file__, '../..'))\n",
            "pytpath = os.path.join(tree, 'PyTransport', 'PyTransport')\n",
            "sys.path.append(pytpath)\n",
            "import PyTransSetup\n",
            "nF = {}\n".format(N),
            "nP = {}\n".format(N),
            "f  = sym.symarray('f',nF)\n",
            "p  = sym.symarray('p',nP)\n",
            "s  = [sym.Rational(1,2) * f[i] ** 2 * p[i] ** 2 for i in range({N})]\n".format(N=N),
            "V  = sum(s)\n",
            "PyTransSetup.tol(1E-8,1E-8)\n"
            "PyTransSetup.potential(V,nF,nP,silent=False)\n",
            "PyTransSetup.compileName('{}Quad')".format(N)]


# Finally we write lines to file
with open(fname, "w") as f:
    for line in preamble:
        f.writelines(line)

