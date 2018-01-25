#!/usr/bin/env python
import os
from os.path import split, join, realpath
import sys
import numpy as np
import shutil
#import time
import subprocess

#import computevoltage_ForHLR as cv
import computeVoltage_massProd as cv


''' call that script via: python example.py *json path_to_TMP
    It will read in the information in the jso-file, event-by-event- and convert it to the nessary format of input for the radio morphing.
    For every event it creates a subfolder in $TMP/InterpolatedSignals/ having the name set by the event-tag for a later identification of the event. 
    The script produces teh needed antpos.dat in the event-folder, caclulated traces after interpolation saved in the same folder
''' 

## VOLTAGE COMPUTAION: if need set to 1
VOLTAGE=1


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


json_file=str(sys.argv[1])
#print "json file : ", json_file

### move json file to tmp_dir
# get name and path
path_json, filename = os.path.split(json_file)
#mv json file to TMP
shutil.copy(json_file,tmp_dir)
# get new path and json file 
json_file = join(tmp_dir, filename) # original json file containing a bunch of events
#print json_file


#j=0
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
        
        
        #### Neutrino energy - of the regenerated neutrino, not the primary one
        #print "nu momentum: ", event["decay"][0][1]
        #nu_energy=np.sqrt((event["decay"][0][1][0])**2 + (event["decay"][0][1][1])**2 + (event["decay"][0][1][2])**2)* 1e-9 # GeV in EeV
               

        ####### STUDY IMPACT OF SEVERAL DECAY PRODUCTS 
        subshowers=0. # subshowers introduced for each decayproduct, subshower=0. leading particle gets all the energy
        
        if subshowers==0:
                ### ENERGY
                ep_array=np.zeros(len(event["decay"])-1)
                for i in range(1, len(event["decay"])): #len(event["decay"])-1 # gives you the number of decay products in event, i==0: neutrino
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

                # Perform RADIO MORPHING
                radiomorphing.process(sim_dir, shower, antennas, out_dir)
                
                
                #print nu_energy, int(nu_energy), int(nu_energy)*1e18, "{:1.0e}".format(int(nu_energy)*1e18), int(nu_energy*1e18),"{:1.0e}".format(int(nu_energy*1e18))
                
                #set up the you folder structure within data_dir like: latitude-longitude/showerenergy/theta/phi
                # ATTENTION: not clear if nu_energy from momentum correct
                structure=join("X"+str(int(decay_pos[0]))+"_Y"+str(int(decay_pos[1])), "E{:1.0e}".format(int(ep*1e18)), "T"+str(int(theta)), "P"+str(int(azimuth)) ) #"{:1.0e}".format(int(nu_energy*1e18))
                print structure
                structure=join(data_dir, structure)
                if not os.path.exists(structure): # later this is not necessary with $TMP
                    os.makedirs(structure)
                        

                
                ##### VOLTAGE COMPUTATION
                if VOLTAGE==1:
                    #### get VOLTAGE traces out_#.txt in out_dir
                    #alpha_sim=0 # ATTENTIONhas to be handed over as an array at some point
                    effective = 1 # use effective zenith
                    
                    #cv.compute(out_dir, alpha_sim, effective, json_file)
                    #### cv produces a new jsonfile named eventtag.voltage.json in tmp_dir for every event in json_file
                    cv.compute('json',out_dir, effective,theta, azimuth, ep, height, primary,json_file)
                    
                    cvjson_file=join(tmp_dir, "InterpolatedSignals",str(event["tag"])+".voltage.json") # name of new created jsonfile for each event, folder level same as for event folder
                    
                    
                    #newname= join(data_dir, str(event["tag"])+".voltage.json")
                    #shutil.copy(cvjson_file,newname)
                    shutil.move(cvjson_file, structure)
                    print "Move json file "+     str(event["tag"])+".voltage.json"       +" moved to ",    structure
                
                
                
                import tarfile
                tar_name= join(tmp_dir, "InterpolatedSignals", str(event["tag"])+".tgz")
                tar = tarfile.open(tar_name, "w:gz")
                tar.add(out_dir, arcname=str(event["tag"]))
                tar.close()
                
                
                # copy out_dir from $TEMP to $PROJECT (data_dir), rm out_dir
                #shutil.move(out_dir, data_dir) 
                print tar_name
                shutil.move(tar_name, structure) 
                shutil.rmtree(out_dir)
                
                
                
                
                ### Upload to iRODS: write a shell script and start it
                #import subprocess
                #subprocess.call(['./test.sh'])

                tgzfile=structure+"/"+str(event["tag"])+".tgz"
                jfile=structure+"/"+str(event["tag"])+".voltage.json"

                sh_file= tmp_dir + "upload.ish"
                print sh_file
                inpfile = open(sh_file,"w+")
                
                inpfile.write("#!/project/fh1-project-huepra/le6232/soft/ishell/bin/ishell")
                inpfile.write("# Example script for copying data from fh1 to iRODS. Note")
                inpfile.write("# that you need to set your iRODS environment first.")
                inpfile.write("#")
                inpfile.write("# -----------------------------------------")
                inpfile.write("# Scheduler options")
                inpfile.write("# -----------------------------------------")
                inpfile.write("#MSUB -N put-events")
                inpfile.write("#MSUB -q singlenode")
                inpfile.write("#MSUB -l nodes=1:ppn=1")
                inpfile.write("#MSUB -l walltime=24:00:00")
                inpfile.write("#MSUB -l pmem=3gb")
                inpfile.write("#MSUB -v PATH, PYTHONPATH")
                inpfile.write("# -----------------------------------------")
                inpfile.write("put -f {:s}\n".format(tgzfile))
                inpfile.write("put -f {:s}\n".format(jfile))
                inpfile.write("#")
                
                inpfile.close()
                
                subprocess.call([sh_file])
                
                os.remove(sh_file)
                
                
                
                
                
print "Job done"





#### NOTE Running subshowers commented out the the moment since timing not yet included. Has to be updated if evertyhing is included
        #if subshowers==1:
            ##### here if necessary loop over the single decay products. Create a folder for ecah of them with the event-tag_i for a later summing up of the traces.
            ##### check whether running subshower or leading product with total energy (like done in ToyModel) necessary
            ##### Then, move it after antenna and take RM into the loop
            
            #for i in range(1, len(event["decay"])): #len(event["decay"])-1 # gives you the number of decay products in event
                #if float(event["decay"][i][0]) in particle_list: # just for valid particles
                    #pp=event["decay"][i][1] # momentum vector, second decay product: event["decay"][2][1] 
                    #ep=np.sqrt(pp[0]**2+pp[1]**2+pp[2]**2)# in GeV
                    #ep*=1.e-9 #in EeV
                    #print "energy in EeV: ", ep
                    
                    #### PID primary
                    ##part_dic={'221.0':'eta','211.0': 'pi+', '-211.0': 'pi-','111.0': 'pi0', '22.0':'gamma', '13.0':'muon', '11.0': 'electron', '15.0':'tau', '16.0':'nu(t)', '321.0': 'K+', '-321.0': 'K-','130.0':'K0L', '310.0':'K0S','-323.0':'K*+'}
                    #PID= float(event["decay"][i][0])
                    #el_list=[22.0,11.0,-11.0, 111.0] #'22.0':'gamma', '11.0': 'electron', '-11.0' positron: '111.0': 'pi0'
                    #if PID in el_list:
                        #primary="electron"
                    #else: # pion-like
                        #primary="pion"
                    #print primary
                    
                    
                    ### create a folder in $TMP for each event
                    #out_dir = join(tmp_dir, "InterpolatedSignals", str(event["tag"])) # will be created in $TMP, will be deleted after each event
                    #out_dir=out_dir+"_"+str(i) # folder for each decay particle
                    #if not os.path.exists(out_dir):
                        #os.makedirs(out_dir)
                    #print "folder ", out_dir 
                    
                
                ###save ANTENNA POSITIONS in anpos.dat for each event
                ### antenna positions have to be corrected for the deacy positions since decay in morphing at (0,0,height) 
                    #antennas = join(out_dir, "antpos.dat") # file of list of desired antennas position -> in out_dir at $TMP
                    #correction= np.array([decay_pos[0], decay_pos[1], 0.])
                    #print "correction: ",correction
                    #print "antenna: ", event["antennas"][0]
                    #np.savetxt(antennas, event["antennas"]-correction, delimiter='  ',fmt='%.1f')   # in GPS coordinates, later being read in in core.py
                    #print "antenna corrected: ", event["antennas"][0]-correction
                    

                    ###### Start radio morphing

                    #shower = {
                            #"primary" : primary,       # reference shower "electron" at the moment
                            #"energy" : ep,               # EeV
                            #"zenith" : theta,               # deg (GRAND frame)
                            #"azimuth" : azimuth,                # deg (GRAND frame)
                            #"injection_height" : height ,    # m
                            #"altitude" : decay_altitude}   # m

                        ## Perform the radiomorphing
                    #radiomorphing.process(sim_dir, shower, antennas, out_dir)
                    
                    ##NOTE: traces of shubshowers have to be added up for comparison
                    
                    
                    ### copy out_dir from $TEMP to $PROJECT (data_dir), rm out_dir
                    ##shutil.move(out_dir, data_dir) 
                    
                    #import tarfile
                    #tar_name= join(tmp_dir, "InterpolatedSignals", str(event["tag"])+".tgz")
                    #tar = tarfile.open(tar_name, "w:gz")
                    #tar.add(out_dir, arcname=str(event["tag"]))
                    #tar.close()
                    
                    ##import tarfile

                    ##with tarfile.open( out_dir + ".tgz", "w:gz" ) as tar:
                        ##for name in os.listdir( out_dir):
                            ##tar.add(name)
                            ##print "tar-file: ", name
                    
                    
                    ## copy out_dir from $TEMP to $PROJECT (data_dir), rm out_dir
                    ##shutil.move(out_dir, data_dir) 
                    #print tar_name
                    #shutil.move(tar_name, data_dir) 
                    #shutil.rmtree(out_dir)
        

#if VOLTAGE==1:
    #filename=os.path.splitext(filename)[0]
    #newname= join(data_dir, filename+".voltage.json")
    #shutil.copy(json_file,newname)
    #print "Move json file to:", newname 


#try:
    #shutil.rmtree(tmp_dir) # remove tmp dir, but not necessary on ForHLR, done automatically when jobs finished
    #print "TMP deleted"
#except IOError:
    #print "Job done"

