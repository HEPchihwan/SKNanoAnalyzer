#!/bin/bash

echo "Setting up environment for TTbar plotting..."

# Source the SKNanoAnalyzer setup
cd /data6/Users/achihwan/SKNanoAnalyzer
source setup.sh

# Go back to plotting directory
cd /data6/Users/achihwan/SKNanoAnalyzer/plots/TTBar/2022

echo "Running TTbar Data vs MC plotter..."
python ttbar_data_mc_plotter.py