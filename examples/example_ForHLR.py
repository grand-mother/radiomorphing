#!/usr/bin/env python
import os
from os.path import split, join, realpath
import sys
import numpy as np
import shutil
#import time


''' call that script via: python example.py *json path_to_TMP
    It will read in the information in the jso-file, event-by-event- and convert it to the nessary format of input for the radio morphing.
    For every event it creates a subfolder in $TMP/InterpolatedSignals/ having the name set by the event-tag for a later identification of the event. 
    The script produces teh needed antpos.dat in the event-folder, caclulated traces after interpolation saved in the same folder
''' 

# Expand the PYTHONPATH and import the radiomorphing package #NOTE: this would be on the shared disc
root_dir = realpath(join(split(__file__)[0], "..")) # = $PROJECT
sys.path.append(join(root_dir, "lib", "python"))

import radiomorphing
import retro



EARTH_RADIUS=6370949. #m

## Settings of the radiomorphing
data_dir = join(root_dir, "examples", "data") # will be in shared storage $PROJECT
simref_dir = join(data_dir, "GrandEventADetailed2") # will be copied from $PROJECT to $TMP, scaled_Evet* can be stored there since it will be overwritten for each event

# ATTENTION: set here the path to the temporary storage at evry core, for testing set to test_dir created
tmp_dir=sys.argv[2] #join(root_dir, "examples", "test_temp") # "$TMP" at core

if not os.path.exists(tmp_dir): # later this is nit necessary with $TMP
    os.makedirs(tmp_dir)
    print "PATH to TMP ", tmp_dir 
sim_dir=join(tmp_dir, "GrandEventADetailed2") # local copy of refernce shower at $TMP, plus folder which will contain the scaled traces 


#t0=time.time()
if not os.path.exists(sim_dir):
    shutil.copytree(simref_dir, sim_dir)#tmp_dir+"/GrandEventADetailed2") # copy refernce shower to TMP of core
print "path to sim " , sim_dir
#t1=time.time()
#print " time needed for copy :", t1-t0
    
particle_list=[22.0, 11.0, -11.0, 111.0, 211.0, -211.0, 221.0] # 22:gamma, 11:e+-, 111:pi0, 211:pi+-, 211:eta

print "json file : ", str(sys.argv[1])
#j=0
### MAYBE: this has to be done in a script which is one level higher and calling the example.py
from retro.event import EventIterator
for event in EventIterator(str(sys.argv[1])):#"events-flat.json"): #json files contains a list of events which shall run on one node"
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
   #j=j+1
   #if j<43: 
        print "\n"
        print "Event ", str(event["tag"]), " started"
                            
                            
        ###DECAY
        decay_pos=event["tau_at_decay"][1]
        height=decay_pos[2]
        print "decay position: ", decay_pos
        decay_pos=decay_pos+np.array([0.,0.,EARTH_RADIUS]) # corrected for earth radius
        print "decay position after correction: ", decay_pos
        
        
        decay_altitude=event["tau_at_decay"][3] 
        print "decay decay_altitude: ", decay_altitude
        
        v=event["tau_at_decay"][2]# shower direction, assuming decay products strongly forward beamed  

        
        ###ANGLES
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
        subshowers=0. # subshowers introduced for each decayproduct, subshower=0. leading particle gets all the energy
        
        if subshowers==0:
                ### ENERGY
                ep_array=np.zeros(len(event["decay"])-1)
                for i in range(1, len(event["decay"])): #len(event["decay"])-1 # gives you the number of decay products in event
                    if float(event["decay"][i][0]) in particle_list: # just for valid particles 
                        pp=event["decay"][i][1] # momentum vector, second decay product: event["decay"][2][1] 
                        ep_array[i-1]=np.sqrt(pp[0]**2+pp[1]**2+pp[2]**2)# in GeV
                    print "particle ", str(i), "PID:",event["decay"][i][0]," energy in EeV: ", ep_array[i-1]*1e-9 #np.sqrt(pp[0]**2+pp[1]**2+pp[2]**2)* 1.e-9 
                ep= np.sum(ep_array)* 1.e-9 # GeV in EeV
                print "energy in EeV: ", ep
                
                ### PID primary
                #part_dic={'221.0':'eta','211.0': 'pi+', '-211.0': 'pi-','111.0': 'pi0', '22.0':'gamma', '13.0':'muon', '11.0': 'electron', '15.0':'tau', '16.0':'nu(t)', '321.0': 'K+', '-321.0': 'K-','130.0':'K0L', '310.0':'K0S','-323.0':'K*+'}
                particle=int(np.argmax(ep_array) +1) # not forget about the inital neutrino and array start with 0
                PID= float(event["decay"][int(np.argmax(ep_array) +1)][0])
                el_list=[22.0, 11.0, -11.0, 111.0] #'22.0':'gamma', '11.0': 'electron', '-11':positron,  '111.0': 'pi0'
                if PID in el_list:
                    primary="electron"
                else: # pion-like
                    primary="pion"
                print primary
                
                
                ## create a folder in $TMP for each event
                out_dir = join(tmp_dir, "InterpolatedSignals", str(event["tag"])) # will be created in $TMP, will be deleted after each event 
                if not os.path.exists(out_dir):
                    os.makedirs(out_dir)
                print "folder ", out_dir 
                
            
            ##save ANTENNA POSITIONS in anpos.dat for each event
                antennas = join(out_dir, "antpos.dat") # file of list of desired antennas position -> in out_dir at $TMP
                correction= np.array([decay_pos[0], decay_pos[1], 0.])
                print "correction: ",correction
                print "antenna: ", event["antennas"][0]
                np.savetxt(antennas, event["antennas"]-correction, delimiter='  ',fmt='%.1f')   # in GPS coordinates
                print "antenna corrected: ", event["antennas"][0]-correction
                

                ##### Start radio morphing

                shower = {
                        "primary" : primary,       # reference shower "electron" at the moment
                        "energy" : ep,               # EeV
                        "zenith" : theta,               # deg (GRAND frame)
                        "azimuth" : azimuth,                # deg (GRAND frame)
                        "injection_height" : height,    # m
                        "altitude" : decay_altitude}   # m

                    # Perform the radiomorphing
                radiomorphing.process(sim_dir, shower, antennas, out_dir)
                
                
                import tarfile
                tar_name= join(tmp_dir, "InterpolatedSignals", str(event["tag"])+".tgz")
                tar = tarfile.open(tar_name, "w:gz")
                tar.add(out_dir, arcname=str(event["tag"]))
                tar.close()
                
                #import tarfile

                #with tarfile.open( out_dir + ".tgz", "w:gz" ) as tar:
                    #for name in os.listdir( out_dir):
                        #tar.add(name)
                        #print "tar-file: ", name
                
                
                # copy out_dir from $TEMP to $PROJECT (data_dir), rm out_dir
                #shutil.move(out_dir, data_dir) 
                print tar_name
                shutil.move(tar_name, data_dir) 
                shutil.rmtree(out_dir)


        if subshowers==1:
            #### here if necessary loop over the single decay products. Create a folder for ecah of them with the event-tag_i for a later summing up of the traces.
            #### check whether running subshower or leading product with total energy (like done in ToyModel) necessary
            #### Then, move it after antenna and take RM into the loop
            
            for i in range(1, len(event["decay"])): #len(event["decay"])-1 # gives you the number of decay products in event
                if float(event["decay"][i][0]) in particle_list: # just for valid particles
                    pp=event["decay"][i][1] # momentum vector, second decay product: event["decay"][2][1] 
                    ep=np.sqrt(pp[0]**2+pp[1]**2+pp[2]**2)# in GeV
                    ep*=1.e-9 #in EeV
                    print "energy in EeV: ", ep
                    
                    ### PID primary
                    #part_dic={'221.0':'eta','211.0': 'pi+', '-211.0': 'pi-','111.0': 'pi0', '22.0':'gamma', '13.0':'muon', '11.0': 'electron', '15.0':'tau', '16.0':'nu(t)', '321.0': 'K+', '-321.0': 'K-','130.0':'K0L', '310.0':'K0S','-323.0':'K*+'}
                    PID= float(event["decay"][i][0])
                    el_list=[22.0,11.0,-11.0, 111.0] #'22.0':'gamma', '11.0': 'electron', '-11.0' positron: '111.0': 'pi0'
                    if PID in el_list:
                        primary="electron"
                    else: # pion-like
                        primary="pion"
                    print primary
                    
                    
                    ## create a folder in $TMP for each event
                    out_dir = join(tmp_dir, "InterpolatedSignals", str(event["tag"])) # will be created in $TMP, will be deleted after each event
                    out_dir=out_dir+"_"+str(i) # folder for each decay particle
                    if not os.path.exists(out_dir):
                        os.makedirs(out_dir)
                    print "folder ", out_dir 
                    
                
                ##save ANTENNA POSITIONS in anpos.dat for each event
                ## antenna positions have to be corrected for the deacy positions since decay in morphing at (0,0,height) 
                    antennas = join(out_dir, "antpos.dat") # file of list of desired antennas position -> in out_dir at $TMP
                    correction= np.array([decay_pos[0], decay_pos[1], 0.])
                    print "correction: ",correction
                    print "antenna: ", event["antennas"][0]
                    np.savetxt(antennas, event["antennas"]-correction, delimiter='  ',fmt='%.1f')   # in GPS coordinates, later being read in in core.py
                    print "antenna corrected: ", event["antennas"][0]-correction
                    

                    ##### Start radio morphing

                    shower = {
                            "primary" : primary,       # reference shower "electron" at the moment
                            "energy" : ep,               # EeV
                            "zenith" : theta,               # deg (GRAND frame)
                            "azimuth" : azimuth,                # deg (GRAND frame)
                            "injection_height" : height ,    # m
                            "altitude" : decay_altitude}   # m

                        # Perform the radiomorphing
                    radiomorphing.process(sim_dir, shower, antennas, out_dir)
                    
                    #NOTE: traces of shubshowers have to be added up for comparison
                    
                    
                    # copy out_dir from $TEMP to $PROJECT (data_dir), rm out_dir
                    shutil.move(out_dir, data_dir) 
        

#try:
    #shutil.rmtree(tmp_dir) # remove tmp dir, but not necessary on ForHLR, done automatically when jobs finished
    #print "TMP deleted"
#except IOError:
print "Job done"

