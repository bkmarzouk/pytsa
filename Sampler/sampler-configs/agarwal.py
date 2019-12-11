"""

To build a sampler, we are required to write a configuration file as follows

"""

# Standard preamble of modules
import numpy as np
import os
import sys


""" Computation Configuration | CHANGE ME """


sampler_name = "simple_agarwal_new"
save_location = "default" # If default: outputs will be written to sampler-build directory


sampler = {
    "SamplerName": sampler_name,          # Sampler name, should be unique from other sampler installations
    "PyTransportModule": "PyTransagarwal_dmax_6pt0",  # Installed module, e.g. PyTrans2Quad
    "NSamples": 1,                      # Compute N samples that support minimum number of efolds (see accept crit.)
    "ExitTime": 55,                       # Assign number of efoldings before end of inflation that scales exit horizon
    "SubEvolution": 6,                    # Duration (efoldings) of sub-horizon evolution
    "tols": [1e-8, 1e-8]                  # Integration tolerance
}

""" Define system paths on construction | DO NOT CHANGE """
system = {
    # sys01
    # sys02
    # sys03
}

computations = {
    "2pf": False,  # Compute 2-point function observables
    "3pf": False,  # Compute 3-point function observables
    "Mij": False   # Evaluate Mass-Matrix Eigenvalues @ ExitTime, (indirectly) required if testing adiabatic limit
}


which3pf = [
    {"config_name": "eq", "latex": "f_{NL}^{eq}", "alpha": 1./3., "beta": 1./3.}
]


accept_criteria = { # Accepts trajectories subject to...
    "MinimumEfolds": 60,       # Required: Minimum amount of inflation for successful sampler
    "TestAdiabaticLimit": False # Not Required: Infers adiabatic limit from background evolution
}


end_inflation = { # End inflation under slow-roll violation. If false, uses field conditions to end realisation
    "Canonical": False,    # \epsilon > 1
}


end_conditions = [ # If non-canonical end, specify field relevant field values
    {
        "FieldNumber": 0,
        "FieldValue" : 0.02,
        "FromAbove"  : True,
        "Successful" : True
    },

    {
        "FieldNumber": 0,
        "FieldValue" : 1.0,
        "FromAbove"  : False,
        "Successful" : False
    }

]



""" Set Sampling Strategy, i.e. priors on initial conditions and parameters | CHANGE ME """



required_modules = [ # Specify modules required for sampling strategy
    "numpy as np"
]

field_positions = [ # Set initial field-space positions
    {
        "FieldNumber": "ALL",
        "Command": 1.0
    },
    {
        "FieldNumber": 0,
        "Command": 0.9
    }
]

field_velocities = [
    {
        "FieldNumber": "ALL",
        "Command": 0
    }
]

parameter_values = [ # If None is given for latex, parameter will *not* be logged
    {
        "ParameterNumber": "ALL",
        "Command": "np.random.normal(0, 1)",
        "LaTeX": None
    },
    
    {
        "ParameterNumber": 0,
        "Command": "np.random.normal(1e-2, 1e-3)",
        "LaTeX": "T_3"
    },
    
    {
        "ParameterNumber": 1,
        "Command": "np.random.normal(1e-3, 1e-4)",
        "LaTeX": "a_0"
    },
    
    {
        "ParameterNumber": 2,
        "Command": "np.random.uniform(0, 1)",
        "LaTeX": "Q"
    },
    
    {
        "ParameterNumber": 3,
        "Command"        : "np.random.normal(1e-1, 1e-2)",
        "LaTeX"          : "\phi_{uv}"
    },
    
    {
        "ParameterNumber": 4,
        "Command"        : "np.random.normal(0, 1e-14)",
        "LaTeX"          : "V_0"
    },
]
