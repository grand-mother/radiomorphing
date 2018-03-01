# Radiomorphing
Welcome to Radio Morphing!

These people made that tool amazing:
W. Carvalho, K. de Vries, O. Martineau, C. Medina, V. Niess, M. Tueros, A. Zilles

Details of the methods can be found in <RM_PUBLICATION>

===== This version is the very first draft and more like a list of notes ====

===== beta users should already find all needed informations to run RM ====


## Description/Introduction

Preliminary example of a radiomorphing Python package.

This python package allows you to calculates the radio signal of an air shower with desired shower parameters which will be detected by an radio array within a fraction of the time common air-simulation tools will need.

The parametrisation was developed for 


## Installation

Coming soon ... E.g. it could be packaged with the standard setup.py tools.

### Run the example
To run the example, execute 
`python examples/example.py`

Following information have to be handed over in that file:
The desired shower parameter as the primary type (electron or pion), its primary energy, its injection height and the injection altitude. Both are needed since for very inclined showers and depending on the definition of the coordinate origin they must not be equivalent due to Earth's curvature. 
In addition, one has to hand over the direction of propagation (zenith and azimuth of the shower in [GRAND conventions](https://github.com/grand-mother/simulations/blob/master/GRANDAngularConventions.pdf)).

The script will read in the antenna positions and electric field traces of the example reference shower given in `examples/data/GrandEventADetaild2` and write the output electric field traces for the desired antenna positions defined in `antpos_desired.dat` to `examples/data/InterpolatedSignals`

### The example reference shower
The example reference shower is an air shower induced by an electron of an energy of 1EeV and an height of 2000m above sealevel. The propagation direction is given by a zenith of 89.5deg (upward-going) and an azimuth of 0deg.

## Documentation

See the [example.py](examples/example.py) script for an example of usage.

The basis of the calclulation is an simulated reference shower. At the moment just results of ZHAireS simulations can be read in. A usage of CoREAS output will be integrated asap.

The internal coordinate system used in radiomorphing is defined the injection point of the shower which is given by (0,0,height). The altitude wrt the sealevel of the injection has to be handed over as well for the scaling part of the method. 
In comparison to ZHAireS or CoREAS simulations, Radio Morphing expects the direction of propagation in  [GRAND conventions](https://github.com/grand-mother/simulations/blob/master/GRANDAngularConventions.pdf) as input.

 -explain units for all parameters and outputs
 -explain scaling (shortly)
 -explain interpolation (shortly) -> doesn't work for very low/high frequencies, define the valid band (upsampling +downsampling done within the script) 
 
 - Magnetic field values are hardcoded for the GRAND example location (Ulastai) in [scaling.py](https://github.com/grand-mother/radiomorphing/blob/afc77779bb0acc09e3458e9e5f0c0e68b77456f9/lib/python/radiomorphing/scaling.py#L287-L291). If needed that values have to be adapted.


### Input overview
to be handed over as in [example.py](https://github.com/grand-mother/radiomorphing/blob/master/examples/example.py)

**set of shower parameters** 
-primary type (electron/pion), 
-energy of the particle inducing the shower in EeV, 
-zenith in deg (GRAND conv), 
-azimuth in deg (GRAND conv), 
-injection height in m== z component of the injection point in m, used to define the injectionposition as (0,0,injectionheight) as reference
-altitude == actual height in m of injection above sealevel which can differ from the injectionheight in ase of Earth's curvature and differing original coordinate system for the desired antenna positions

**desired antenna positions**: list of desired antenna positions handed over in x (along North-South), y (along East-West), z (vertical, r.t. sealevel), since positions must be given in the reference frame defined by injection point, for example saved like  [antpos_desired.dat](https://github.com/grand-mother/radiomorphing/blob/master/examples/data/InterpolatedSignals/antpos_desired2.dat)
 
**path reference shower**: the reference shower is simulated in a star shape pattern (see folder GrandEventADetailed2, -> ask me for the 16 planes)


### Output
out_dir: folder for output
-a#.trace: electric field traces for EW,NS, and UP component, time in ns, electric field strength in muV/m (cross-check that), #=ID number of antenna, meaning index of antenna in teh antenna list

### Setting up own study/referenc shower
 -> here there should come the explanation how to set up an own reference shower
 ecah plane needs to be stored in its own subfolder (all traces and corresponding antenna positions in an axtra file)
 


## Future projects
 - 
 - include CoREAS simus as possible input 
 

## License

The radiomorphing package is under the **GNU LGPLv3** license. See the provided
[LICENSE](LICENSE) and [COPYING.LESSER](COPYING.LESSER) files.
