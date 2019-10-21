"""

To build a sampler, we are required to write a configuration file as follows

"""

# Standard preamble of modules
import numpy as np
import os
import sys
import dir_structure as ds



""" Path Configuration | AUTOMATIC """


# Get python path: This should be automatic assuming directory structure has not been altered
pytpath, savedir = ds.get_pyt_paths()


""" Computation Configuration | CHANGE ME """


sampler_name = "dquad_example"
saveloc = os.path.join(savedir, sampler_name)


sampler = {
    "SamplerName": sampler_name,          # Sampler name, should be unique from other sampler installations
    "SaveLocation": saveloc,              # Save location for outputs
    "PyTransportModule": "PyTrans2Quad",  # Installed module, e.g. PyTrans2Quad
    "NSamples": 200,                      # Compute N samples that support minimum number of efolds (see accept crit.)
    "ExitTime": 55,                       # Assign number of efoldings before end of inflation that scales exit horizon
    "SubEvolution": 6,                    # Duration (efoldings) of sub-horizon evolution
    "tols": [1e-8, 1e-8]                  # Integration tolerance
}


computations = {
    "2pf": True,  # Compute 2-point function observables
    "3pf": True,  # Compute 3-point function observables
    "Mij": True   # Evaluate Mass-Matrix Eigenvalues @ ExitTime, (indirectly) required if testing adiabatic limit
}


which3pf = [
    {"config_name": "eq", "latex": "f_{NL}^{eq}", "alpha": 1./3., "beta": 1./3.}
]


accept_criteria = { # Accepts trajectories subject to...
    "MinimumEfolds": 60,       # Required: Minimum amount of inflation for successful sampler
    "TestAdiabaticLimit": True # Not Required: Infers adiabatic limit from background evolution
}


end_inflation = { # End inflation under slow-roll violation. If false, uses field conditions to end realisation
    "Canonical": True,    # \epsilon > 1
}


end_conditions = [ # If non-canonical end, specify field relevant field values
    None
]



""" Set Sampling Strategy, i.e. priors on initial conditions and parameters | CHANGE ME """



required_modules = [ # Specify modules required for sampling strategy
    "numpy", "scipy"
]

field_positions = [ # Set initial field-space positions
    {"FieldNumber": "ALL",
     "Command": "numpy.random.normal(-20, 20)"
    },
    {"FieldNumber": 0,
     "Command": "numpy.random.uniform(-10, 10)"
    }
]

field_velocities = [
    {"FieldNumber": "ALL",
     "Command": "SlowRoll"
    }
]

parameter_values = [ # If None is given for latex, parameter will *not* be logged
    {"ParameterNumber": 0,
     "Command": 1e-3,
     "LaTeX": "m_\chi"
    },

    {"ParameterNumber": 1,
     "Command": "numpy.exp(numpy.random.uniform(-2, -5))",
     "LaTeX": "m_\phi"
    }
]
