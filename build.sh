#!/usr/bin/env bash

set -e  
set -u  

PYTHON=${PYTHON:-python3}       
VENV_DIR=".venv"                
REQ_FILE="requirements.txt"         

echo ">>> Building Python project..."
echo ">>> Using Python: $PYTHON"
echo ">>> Virtual environment directory: $VENV_DIR"

if [ ! -d "$VENV_DIR" ]; then
    echo ">>> Creating virtual environment..."
    $PYTHON -m venv "$VENV_DIR"
else
    echo ">>> Virtual environment already exists."
fi

source "$VENV_DIR/bin/activate"

echo ">>> Upgrading pip..."
pip install --upgrade pip setuptools wheel

if [ -f "$REQ_FILE" ]; then
    echo ">>> Installing dependencies from $REQ_FILE..."
    pip install -r "$REQ_FILE"
else
    echo "!!! WARNING: $REQ_FILE not found. Skipping dependency installation."
fi

echo ">>> Build complete!"
echo "To activate the virtual environment later, run:"
echo "    source $VENV_DIR/bin/activate"
