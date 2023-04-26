#!/bin/bash

# Install packages with pip
python -m venv dependencies
source dependencies/bin/activate
pip install -r requirements.txt

# Install packages with apt-get using sudo
#sudo apt-get update
#sudo apt-get install -y GDAL

# Add any other installation commands as needed

chmod +x Eagloo.desktop
chmod +x test_launch.sh
