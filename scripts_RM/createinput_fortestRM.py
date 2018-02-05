#!/usr/bin/env python


''' call that script via: python createinput_fortestRM *json 
    It will read in the information in the json-file, event-by-event- and convert it to the nessary format of input for starting Zhaires sims to cross-check output of RM.
    One can distinguish/test subshowers or leading-particle showers by setting subshower=0/1
''' 



import os
from os.path import split, join, realpath
import sys
import numpy as np
import shutil
#import time
import random


# Expand the PYTHONPATH and import the radiomorphing package #NOTE: this would be on the shared disc
root_dir = realpath(join(split(__file__)[0], "..")) # = $PROJECT
sys.path.append(join(root_dir, "lib", "python"))

import retro

EARTH_RADIUS=6370949. #m
GEOMAGNET = (56.5, 63.18, 2.72) # Geomagnetic field (Amplitude [uT], inclination [deg], declination [deg]).



'''NOTE: to be set to 1 (=create subshower) or 0 (=shower by leading particle)'''
subshowers=0. # subshowers introduced for each decayproduct, subshower=0. leading particle gets all the energy




##########################################################################################################
def generate_input(task=0,energy=None, azimuth=None, zenith=None, products=None, height=None, antennas=None, info=None):
    """Generate the input stream for ZHAireS."""

    zen,azim = GRANDtoZHAireS(zenith,azimuth)

    a=" ".join(map(str, products))
    b="".join( c for c in a if  c not in "(),[]''")

    seed = random.uniform(0.,1.)

    # Format the stream.
    stream = [
        "# {:s}".format(info),
        "AddSpecialParticle      RASPASSProton    /home/renault/zhaires/RASPASSprimary/RASPASSprimary Proton",
        "AddSpecialParticle      RASPASSIron      /home/renault/zhaires/RASPASSprimary/RASPASSprimary Iron",
        "AddSpecialParticle      RASPASSelectron  /home/renault/zhaires/RASPASSprimary/RASPASSprimary Electron",
        "AddSpecialParticle      RASPASSpi+  /home/renault/zhaires/RASPASSprimary/RASPASSprimary pi+",
        "AddSpecialParticle      RASPASSpi-  /home/renault/zhaires/RASPASSprimary/RASPASSprimary pi-",
        "AddSpecialParticle      RASPASSpi0  /home/renault/zhaires/RASPASSprimary/RASPASSprimary pi0",
        "AddSpecialParticle      RASPASSMulti /home/renault/zhaires/RASPASSprimary/RASPASSprimary {:s}".format(b),
        "#########################",
        "TaskName {:s}".format(task),
        "PrimaryParticle RASPASSMulti",
        "PrimaryEnergy {:.5E} EeV".format(energy),
        "PrimaryZenAngle {:.5f} deg".format(zen),
        "PrimaryAzimAngle {:.5f} deg Magnetic".format(azim),
        "ForceModelName SIBYLL",
        "SetGlobal RASPASSHeight {:.5f} m".format(height),
        "RandomSeed {:.5f}".format(seed),
        "########################",
        "PropagatePrimary On",
        "SetGlobal RASPASSTimeShift 0.0",
        "SetGlobal RASPASSDistance 0.00"
    ]

    for a in antennas:    
        stream.append("AddAntenna {:1.2f} {:1.2f} {:1.2f}".format(a[0],a[1],a[2]))

    stream += [
        "##########################",
        "TotalShowers 1",
        "RunsPerProcess Infinite",
        "ShowersPerRun 1",
        "Atmosphere 1",
        "AddSite Ulastai 42.55 deg 86.68 deg {:.3f} m".format(1500.),
        "Site Ulastai",
        "Date 1985 10 26",
        "GeomagneticField On",
        "GeomagneticField {:.4f} uT {:.2f} deg {:.2f} deg".format(*GEOMAGNET),
        "GroundAltitude {:.3f} m".format(0.),
        "ObservingLevels 510 50 g/cm2   900 g/cm2",
        "PerShowerData Full",
        "SaveNotInFile lgtpcles All",
        "SaveNotInFile grdpcles All",
        "RLimsFile grdpcles 0.000001 m 10 km",
        "ResamplingRatio 100",
        "#########################",
        "RLimsTables 10 m 10 km",
        "ELimsTables 2 MeV 1 TeV",
        "ExportTables 5501 Opt a",
        "ExportTables 1293 Opt a",
        "ExportTables 1293 Opt as",
        "ExportTable 1205 Opt a",
        "ExportTable 1205 Opt as",
        "ExportTable 1793 Opt a",
        "ExportTable 1793 Opt as",
        "########################",
        "ForceLowEDecay Never",
        "ForceLowEAnnihilation Never",
        "########################",
        "ZHAireS On",
        "FresnelTime On",
        "FresnelFreq Off",
        "TimeDomainBin 1 ns",
        "AntennaTimeMin -100 ns",
        "AntennaTimeMax 500 ns", #can be extended until 3e-6s but if then there is still nothing then there must be a problem somewhere
        "######################", 
        "ElectronCutEnergy 1 MeV",
        "ElectronRoughCut 1 MeV",
        "GammaCutEnergy 1 MeV",
        "GammaRoughCut 1 MeV",
        "ThinningEnergy 1.e-4 Relative", #It can be 1e-5, 1e-6 or below. But running time inversely proportional to it.
        "ThinningWFactor 0.06"
    ]
    

    return "\n".join(stream)

##########################################################################################################
def GRANDtoZHAireS(zen_DANTON=None, azim_DANTON=0):
    """ Convert coordinates from DANTON convention to ZHAireS convention """

    zen = 180. - zen_DANTON
    azim = azim_DANTON - 180.
    if azim>360:
        azim = azim-360.
    elif azim<0.:
        azim = azim+360.
    return [zen,azim]

##########################################################################################################





    
particle_list=[22.0, 11.0, -11.0, 111.0, 211.0, -211.0, 221.0] # 22:gamma, 11:e+-, 111:pi0, 211:pi+-, 211:eta
part_dic={'221.0':'eta','211.0': 'pi+', '-211.0': 'pi-','111.0': 'pi0', '22.0':'gamma', '13.0':'muon', '11.0': 'electron', '15.0':'tau', '16.0':'nu(t)', '321.0': 'K+', '-321.0': 'K-','130.0':'K0L', '310.0':'K0S','-323.0':'K*+'}

json_file=str(sys.argv[1])
#print "json file : ", json_file



j=0
### MAYBE: this has to be done in a script which is one level higher and calling the example.py
from retro.event import EventIterator
for event in EventIterator(json_file):#"events-flat.json"): #json files contains a list of events which shall run on one node"
   #print event["tag"] # gives you the tag one after eachother, long
   #print event["decay"] #(pid, (momentum_x,, momentum_y, momentum_z)
   #print event["decay"][0] #returns the inital neutrino
   #print event["decay"][1][1][0] # 1st decay product and so on
   #print len(event["decay"])-1 # gives you the number of decay products in event
   #print event["antennas"][0] # greps the first antenna of each event
   #print event["antennas"] # greps all antennas of one event
   #print event["tau_at_decay"] # all information about vertex of decay: energy before decay in GeV, position in the local [x, y, z], direction of propagation, [ux, uy, uz]
   #print event["tau_at_decay"][2]
   
   #### to choose one specific event from a json file or test running on cluster
   j=j+1
   if j<10: 
    #if event["antennas"][0][2]>0:
        print "\n"
        print "Event ", str(event["tag"]), " started"
                            
                            
        ###DECAY
        decay_pos=event["tau_at_decay"][2]
        print "decay position ", decay_pos
        height=decay_pos[2]
        print "decay position: ", decay_pos
        decay_pos=decay_pos+np.array([0.,0.,EARTH_RADIUS]) # corrected for earth radius
        print "decay position after correction: ", decay_pos
        
        
        decay_altitude=event["tau_at_decay"][4][2] 
        print "decay decay_altitude: ", decay_altitude
        
        v=event["tau_at_decay"][2]# shower direction, assuming decay products strongly forward beamed  

        
        ###ANGLES theta, azimuth in deg (GRAND)
        theta = np.degrees(np.arccos(np.dot(v, decay_pos) / np.linalg.norm(decay_pos))) # zenith in GRAND conv.
        print "theta: ", theta
        #orthogonal projection of v onto flat plane to get the azimuth 
        x=np.array([1.,0.,0.]) #NS
        y=np.array([0.,1.,0.]) #EW
        proj_v= np.dot(v,x)*x + np.dot(v,y)*y
#        print proj_v
        azimuth = np.degrees(np.arccos(np.dot(proj_v, x))) # azimuth in GRAND conv., rt NORTH
        if proj_v[1]<0.: # y component of projection negativ, means azimuth >180deg
            azimuth = 360.-azimuth
        print "azimuth: ", azimuth
        

        ####### STUDY IMPACT OF SEVERAL DECAY PRODUCTS 
        #subshowers=1. # subshowers introduced for each decayproduct, subshower=0. leading particle gets all the energy
        
        if subshowers==0: #leading particle gets all the energy
                prefix="lead"
            
                ### ENERGY ep in EeV
                ep_array=np.zeros(len(event["decay"])-1)
                for i in range(1, len(event["decay"])): #len(event["decay"])-1 # gives you the number of decay products in event
                    
                    if float(event["decay"][i][0]) in particle_list: # just for valid particles 
                        pp=event["decay"][i][1] # momentum vector, second decay product: event["decay"][2][1] 
                        ep_array[i-1]=np.sqrt(pp[0]**2+pp[1]**2+pp[2]**2)# in GeV
                    print "particle ", str(i), "PID:",event["decay"][i][0]," energy in EeV: ", ep_array[i-1]*1e-9 #np.sqrt(pp[0]**2+pp[1]**2+pp[2]**2)* 1.e-9 
                etot= np.sum(ep_array)* 1.e-9 # GeV in EeV
                print "energy in EeV: ", etot
                
                ### PID primary - leading particle
                particle=int(np.argmax(ep_array) +1) # not forget about the inital neutrino and array start with 0, index of the primary
                PID= float(event["decay"][int(np.argmax(ep_array) +1)][0]) # the number of particle 
                el_list=[22.0, 11.0, -11.0, 111.0] #'22.0':'gamma', '11.0': 'electron', '-11':positron,  '111.0': 'pi0' - particle has to be in that list
                #if PID in el_list:
                    #primary="electron"
                #else: # pion-like
                    #primary="pion"
                #Get Zhaires name for that particle instead of number code
                #primary ='pi+'#str('pi+')
                #primary =str(PID)
                primary=list(part_dic.values())[list(part_dic.keys()).index(str(PID))]
                print primary
                
                multip=[]
                multip.append((primary, 1.0))
                
                
                print multip

            




        if subshowers==1: # all decay product get simulated
            prefix="sub"
            
            ep_array=[] #np.zeros(len(event["decay"])-1)
            PID_array=[] #np.zeros(len(event["decay"])-1)

            for i in range(1, len(event["decay"])): #len(event["decay"])-1 # gives you the number of decay products in event
                print i, float(event["decay"][i][0])
                if float(event["decay"][i][0]) in particle_list: # just for valid particles
                    pp=event["decay"][i][1] # momentum vector, second decay product: event["decay"][2][1] 
                    ep=np.sqrt(pp[0]**2+pp[1]**2+pp[2]**2)# in GeV
                    ep*=1.e-9 #in EeV
                    print "energy in EeV: ", ep
                    ep_array=np.append(ep_array, ep)
                    
                    ### PID primary
                    #part_dic={'221.0':'eta','211.0': 'pi+', '-211.0': 'pi-','111.0': 'pi0', '22.0':'gamma', '13.0':'muon', '11.0': 'electron', '15.0':'tau', '16.0':'nu(t)', '321.0': 'K+', '-321.0': 'K-','130.0':'K0L', '310.0':'K0S','-323.0':'K*+'}
                    PID= float(event["decay"][i][0])
                    #el_list=[22.0,11.0,-11.0, 111.0] #'22.0':'gamma', '11.0': 'electron', '-11.0' positron: '111.0': 'pi0'
                    #if PID in el_list:
                        #primary="electron"
                    #else: # pion-like
                        #primary="pion"
                    #primary=list(part_dic.values())[list(part_dic.keys()).index(str(PID))]
                    #print primary
                    PID_array=np.append(PID_array,PID)
                    print "PID taken: ", PID
                    
            etot= np.sum(np.asarray(ep_array)) # in EeV  
            print "PID_array ", PID_array, "energy_array ", ep_array
            
            
            multip=[]
            for i in range(0,len(ep_array)):
                multip.append((list(part_dic.values())[list(part_dic.keys()).index(str(PID_array[i]))],ep_array[i]/etot))
            print multip
                    
                    

        ## create a folder in $TMP for each event
        out_dir = "./RMtest_input/"#+ str(event["tag"])
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        print "folder ", out_dir 
                    
                
        ##save ANTENNA POSITIONS in anpos.dat for each event
        ## antenna positions have to be corrected for the deacy positions since decay in morphing at (0,0,height) 
        #antennas = join(out_dir, "antpos.dat") # file of list of desired antennas position -> in out_dir at $TMP
        correction= np.array([decay_pos[0], decay_pos[1], 0.])
        ant=np.copy(event["antennas"])
        ant= np.delete(ant, np.s_[3:5], axis=1)
        #np.savetxt(antennas, event["antennas"]-correction, delimiter='  ',fmt='%.1f')   # in GPS coordinates, later being read in in core.py
        ANTENNAS3=ant-correction
        print "antenna corrected: ", ANTENNAS3
        
        print "attention for flat array correction: z of antenna set to 3m"
        
        # approximation for flat100x100: ant_z=3m (see SLACK with VN)
        for i in range(0,len(ANTENNAS3.T[1])):
                           ANTENNAS3[i,2]=3. #m 
  
  
        path, filename = os.path.split(json_file)
        info= "tag: "+str(event["tag"]) + ", in json:" + filename
  
#### Write the ZHAireS input file
        showerID=str(prefix)+"_"+str(j)
        fileZha=join(out_dir, "SIM_"+str(showerID)+".inp")
        inpfile = open(fileZha,"w+")
        #showerID=str(prefix)+"_"+str(j) # number event in file, tag as comment in .inp-file
        totito  = generate_input(showerID, etot, azimuth, theta, multip, decay_altitude,ANTENNAS3, info)
        inpfile.write(totito)
        inpfile.close()    
        
    #else:
        #print "failed"
                    

print "Job done"


