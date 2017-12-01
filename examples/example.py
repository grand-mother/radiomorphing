#!/usr/bin/env python
import os
import sys

# Expand the PYTHONPATH and import the radiomorphing package
module_path = os.path.join(os.path.split(__file__)[0], "..", "lib", "python")
sys.path.append(os.path.realpath(module_path))
import radiomorphing

# Settings of the radiomorphing
sim_dir = "examples/data/GrandEventADetailed2"
out_dir = "examples/data/InterpolatedSignals"
antennas = "antpos_desired2.dat"

shower = {
    "primary" : "electron",
    "energy" : 0.96,               # EeV
    "zenith" : 89.5,               # deg (GRAND frame)
    "azimuth" : 0.,                # deg (GRAND frame)
    "injection_height" : 2000. }   # m

# Perform the radiomorphing
radiomorphing.process(sim_dir, shower, os.path.join(out_dir, antennas), out_dir)
