#!/bin/bash

CWD="$( cd -- "$( dirname -- "${BASH_SOURCE[0]:-$0}"; )" &> /dev/null && pwd 2> /dev/null; )";

echo $CWD/modelSetup.py

CC=gcc python $CWD/moduleSetup.py install --force

# -- Need to work this out as to not make a mess

#MODELS_DIR=$CWD/../models
#CC=gcc python $CWD/moduleSetup.py example_models --force --prefix $MODELS_DIR