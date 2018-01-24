# Radiomorphing
Welcome to Radio Morphing!

These people made that tool amazing:
W. Carvalho, K. de Vries, O. Martineau, C. Medina, V. Niess, M. Tueros, A. Zilles

Details of the methods can be found in <RM_PUBLICATION>


## Description/Introduction

Preliminary example of a radiomorphing Python package.

This python package allows you to calculates the radio signal of an air shower with desired shower parameters which will be detected by an radio array within a fraction of the time common air-simulation tools will need.



## Installation

Coming soon ... E.g. it could be packaged with the standard setup.py tools.

### Run the example
 How to run the example shower -> include 4 planes of example referenceshower in git and a file containing example positions

## Documentation

See the [example.py](examples/example.py) script for an example of usage.

The basis of the calclulation is an simulated reference shower. 

 -tau decay as reference at 0,0,height - mention reference frame
 -explain angle convention and how they differ from Zhaires and coreas (in case coreas as possible input
 -explain units for all parameters and outputs
 -explain scaling (shortly)
 -explain interpolation (shortly) -> doesn't work for very low/high frequencies, define the valid band (upsampling +donsampling done within the script) 
 


### Input
set of shower parameters: set in example.py
-primary type (electron/pion), 
-energy of the particle inducing the shower in EeV, 
-theta in deg (GRAND conv), 
-azimuth in deg (GRAND conv), 
-injection height == z component of the tau decay in the GRAND study in m, used to define the tau decay position as (0,0,injectionheight)m internally
-altitude == actual height in m of decay above sealevel which can differ from the injectionheight
-antennas: list of desired antenna positions handed over in x (along North-South), y (along East-West), z (vertical, r.t. sealevel), since positions must be given in the reference frame defined by 
 the decay position (where the decay is at (0,0,injectionheight) in m, for example saved as antpos_desired.dat
 
-path reference shower: the refernce shower is simulated in a star shape pattern (see folder GrandEventADetailed2, -> bug me to give you the 16 planes)


### Output
out_dir: folder for output
-a#.trace: electric field traces for EW,NS, and UP component, time in ns, electric field strength in muV/m (cross-check that), #=ID number of antenna, meaning index of antenna in teh antenna list

### Setting up own study
 -> here there should come the explanation how to set up an own reference shower


## Future projects
 - include correct timing so that sub-shower calculation can be peformed (application in methods using the peak time maybe limited)
 - include CoREAS simus as possible input 
 

## License

The radiomorphing package is under the **GNU LGPLv3** license. See the provided
[LICENSE](LICENSE) and [COPYING.LESSER](COPYING.LESSER) files.
