#!/bin/bash

CWD="$( cd -- "$( dirname -- "${BASH_SOURCE[0]:-$0}"; )" &> /dev/null && pwd 2> /dev/null; )";

echo $CWD/modelSetup.py

#CC=gcc python $CWD/moduleSetup.py install --force

MODELS_DIR=$CWD/../models

echo $MODELS_DIR

CC=gcc python $CWD/moduleSetup.py install --force --prefix $MODELS_DIR
