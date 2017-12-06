#!/usr/bin/env python
from os.path import split, join, realpath
import sys

# Expand the PYTHONPATH and import the radiomorphing package
root_dir = realpath(join(split(__file__)[0], ".."))
sys.path.append(join(root_dir, "lib", "python"))
import radiomorphing

# Settings of the radiomorphing
data_dir = join(root_dir, "examples", "data")
sim_dir = join(data_dir, "GrandEventADetailed2")
out_dir = join(data_dir, "InterpolatedSignals")
antennas = join(out_dir, "antpos_desired2.dat")

shower = {
    "primary" : "electron",
    "energy" : 0.96,               # EeV
    "zenith" : 89.5,               # deg (GRAND frame)
    "azimuth" : 0.,                # deg (GRAND frame)
    "injection_height" : 2000. }   # m

# Perform the radiomorphing
radiomorphing.process(sim_dir, shower, antennas, out_dir)
