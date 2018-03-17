
import os
import sys
import numpy as np
from frame import UVWGetter, XYZGetter, get_rotation#, GetUVW
from utils import getCerenkovAngle, load_trace

def _getAirDensity(h):
  #% Using isothermal Model
  rho_0 = 1.225#; % kg/m3
  M = 0.028966#;  %kg/mol
  g = 9.81#; %ms-2
  T = 288.#; % K
  R = 8.32#;
  rho = rho_0*np.exp(-g*M*h/(R*T))
  return rho

###################################

def _getXmax(primarytype, energy, zen2):
    # type of primary (electron or pion, energy in EeV, zen2 (GRAND) in rad
    if primarytype=='electron': # aprroximated by gamma shower
        a=82.5 # g/cm2
        c=342.5 #g/cm2
    if primarytype=='pion': # aprroximated by proton
        a=62.5 # g/cm2
        c=357.5 #g/cm2
    Xmax= a*np.log10(energy*10**6.)+c # E/EeV* 10**6. to be in TeV
    
    #print "for  ", primarytype, " of energy ", energy,  "EeV: Xmax in g/cm2 : ", Xmax

    return Xmax#/abs(np.cos(np.pi-zen2)) # TODO: how to correct for slanted shower

def _dist_decay_Xmax(zen2, injh2, Xmax_primary): #zen2: zenith of target shower
    #% Using isothermal Model
    rho_0 = 1.225*0.001#; % kg/m3 to 0.001g/cm3: 1g/cm3=1000kg/m3, since X given in g/cm2
    M = 0.028966#;  %kg/mol - 1000g/mol
    g = 9.81#; %ms-2
    T = 288.#; % K
    R = 8.32#; J/K/mol , J=kg m2/s2

    hD=injh2
    Xmax_primary= Xmax_primary#* 10. # g/cm2 to kg/m2: 1g/cm2 = 10kg/m2
    gamma=np.pi-zen2 # counterpart of where it goes to
    Re= 6370949 # m, Earth radius
    X=0.
    i=0.
    h=hD
    ai=0
    step=10
    while X< Xmax_primary:
        i=i+1
        ai=i*step #100. #m
        hi= -Re+np.sqrt(Re**2. + ai**2. + hD**2. + 2.*Re*hD - 2*ai*np.cos(gamma) *(Re+hD))## cos(gamma)= + to - at 90dg
        deltah= abs(h-hi) #(h_i-1 - hi)= delta h
        h=hi # new height
        X=X+ rho_0*np.exp(-g*M*hi/(R*T)) * step *100. #(deltah*100) *abs(1./np.cos(np.pi-zen2)) # Xmax in g/cm2, slanted = Xmax, vertical/ cos(theta); density in g/cm3, h: m->100cm, np.pi-zen2 since it is defined as where the showers comes from, abs(cosine) so correct for minus values
    
    #print "Xmax height in m", h, " distance along axis in m ", ai
    return h, ai # Xmax_height in m, Xmax_distance in m

def _scalingfactors(E1, az1, zen1, injh1, E2, az2, zen2, injh2, thetageo, altitude):
    #print "altitude scaling ", altitude
    # 2: target shower, 1: generic shower
    #################### Energy scaling
    #% Energy
    #kE = Esh*1e18./Esh_ref
    kE = E2/E1 # both in 1e18eV

    ############## Azimuth scaling
    #% Azimuth
    Bgeo = [np.sin(thetageo), 0, np.cos(thetageo)]#  % With North = magnetic North
    vref = [np.cos(az1)*np.sin(zen1)*-1., np.sin(az1)*np.sin(zen1)*-1., np.cos(zen1)*-1.]# *-1: account for angle conventions
    vxB_ref = np.cross(vref,Bgeo)
    vxB_ref = np.linalg.norm(vxB_ref)/(np.linalg.norm(vref)*np.linalg.norm(Bgeo))
    v = [np.cos(az2)*np.sin(zen2)*-1., np.sin(az2)*np.sin(zen2)*-1., np.cos(zen2)*-1.]# *-1: account for angle conventions
    vxB = np.cross(v,Bgeo)
    vxB = np.linalg.norm(vxB)/(np.linalg.norm(v)*np.linalg.norm(Bgeo))
    kAz = vxB/vxB_ref   # =cos(az2)/cos(az1)

    #print 'included extra scaling factor cos(Delta azimuth) to respect flipping of components'
    kAz = kAz#* np.cos(az2-az1)

    h_ref=injh1
    h=altitude #injh2 # actual altitude wrt sealevel at decay position of target position 
    #%############## Height+Zenith, distance injection point to xmax rougly 8000m
    #hx_ref = h_ref+8000*np.cos(zen1) #   % Height at reference shower Xmax
    hx_ref = h_ref+8000*np.tan(0.5*np.pi-zen1) #   % Height at reference shower Xmax
    ac_ref = getCerenkovAngle(hx_ref)
    rho_ref = _getAirDensity(hx_ref)
    #hx = h+8000*np.cos(zen2)#   % Height at  shower Xmax
    hx = h+8000*np.tan(0.5*np.pi-zen2)#   % Height at target shower Xmax 
    ac = getCerenkovAngle(hx) # ATTENTION add here altitude to get correct density
    rho = _getAirDensity(hx) # ATTENTION add here altitude to get correct density
    kStretch = float(ac)/float(ac_ref)#  % Stretch factor for the antenna grid
    kRho = np.sqrt(rho_ref/rho)
    kHeight = kRho/kStretch
    kAmp=kE*kAz*kHeight

    #print 'kStretch ', kStretch, ' kAmp ', kAmp,  ' kE ', kE, ' KAz ', kAz, ' kHeight ', kHeight

    return kStretch, kE, kAz, kHeight




def _scalingpulse(dist1, E1, az1, zen1, injh1, E2, az2, zen2, injh2, primary, phigeo, thetageo, l,  positions, path, altitude): # hand over parameters from reference shower 1 and target shower 2, the number of the antenna in the star shape you would like to have  and all position of complete starshape (for strechting needed), the path to the folder containing the sim, and for now the frequencies (should be removed if one included the antenna response

#SCALING
    kStretch, kE, kAz, kHeight = _scalingfactors(E1, az1, zen1, injh1, E2, az2, zen2, injh2, thetageo, altitude)
    kAmp=kE*kAz*kHeight
#    if l==0:
#        print 'kStretch ', kStretch, ' kAmp ', kAmp,  ' kE ', kE, ' KAz ', kAz, ' kHeight ', kHeight

    #print 'Amplitude changed by kAmp= ' + str(kAmp), ' Position changed by kStretch= '+ str(kStretch)



 ###############################################################################
 #### scaling electric fields amplitude
 ################################################

    try:
    ##read in full traces of antenna l: 0:time in ns, 1,2,3: vector field , 4,5,6: efield
        txt1 = load_trace(path, l)[:-2,:]
    except IOError:
        print("antenna ID ",str(int(l)), " file doesn't exist")
        sys.exit()

    
    
    ############## ORIGINAL
    ## Convert efield to shower coordinates to apply the scaling
    #R = get_rotation(zen1, az1, phigeo, thetageo)# original
    #EshowerA = np.dot(txt1[:,1:], R) # original
    
    #### Sciling, kHeight includes 1/kStretch
    #EshowerA.T[0] *= kE * kHeight
    #EshowerA.T[1] *= kE * kAz * kHeight
    #EshowerA.T[2] *= kE * kHeight
    
    ###Backtrafo of efield from shower coord (1,2,3) in xyz (4,5,6) after scaling and/or stretching using the target angles
    #Rt = get_rotation(zen2, az2, phigeo, thetageo).T # original
    #v2 = Rt[:,0]# original
    #txt1[:,1:] = np.dot(EshowerA, Rt)# original
    ############## 
    


######## old school way
    inc=thetageo
    
    az=az1
    zen=zen1


    #same script as for the antenna psoition rotation
    B = np.array([np.cos(phigeo)*np.sin(inc), np.sin(phigeo)*np.sin(inc),np.cos(inc)]) #from oliviers script including phigeo
    B=B/np.linalg.norm(B)
    #print B

    v = np.array([np.cos(az)*np.sin(zen),np.sin(az)*np.sin(zen),np.cos(zen)]) # or *-1: change the direction
    v=v/np.linalg.norm(v)
    #print v

    vxB = np.cross(v,B) #np.array([v[1]*B[2]-v[2]*B[1],v[2]*B[0]-v[0]*B[2],v[0]*B[1]-v[1]*B[0]]) # crossproduct
    vxB = vxB/np.linalg.norm(vxB)
    vxvxB = np.cross(v,vxB) #np.array([v[1]*vxB[2]-v[2]*vxB[1],v[2]*vxB[0]-v[0]*vxB[2],v[0]*vxB[1]-v[1]*vxB[0]])# crossproduct
    vxvxB = vxvxB/np.linalg.norm(vxvxB)
    
    #shape = (len(txt1.T[1]),3)
    EshowerA= np.zeros([len(txt1.T[1]),3])
    
    
    EshowerA.T[0]= txt1.T[1]* v[0] +txt1.T[2]*v[1]+ txt1.T[3]*v[2]
    EshowerA.T[1]= txt1.T[1]* vxB[0] +txt1.T[2]*vxB[1]+ txt1.T[3]*vxB[2]
    EshowerA.T[2]= txt1.T[1]* vxvxB[0] +txt1.T[2]*vxvxB[1]+ txt1.T[3]*vxvxB[2]


    
    ### Scaling, kHeight includes 1/kStretch
    EshowerA.T[0] *= kE * kHeight
    EshowerA.T[1] *= kE * kAz * kHeight
    EshowerA.T[2] *= kE * kHeight

    #print "az2, zen2 ", az2, zen2

    ### angles target shower 
    v2 = np.array([np.cos(az2)*np.sin(zen2),np.sin(az2)*np.sin(zen2),np.cos(zen2)]) # *-1: account for angle conventions
    v2=v2/np.linalg.norm(v2)
    #print v2
    vxB2 = np.cross(v2,B) #np.array([v[1]*B[2]-v[2]*B[1],v[2]*B[0]-v[0]*B[2],v[0]*B[1]-v[1]*B[0]]) # crossproduct
    vxB2 = vxB2/np.linalg.norm(vxB2)
    vxvxB2 = np.cross(v2,vxB2) #np.array([v[1]*vxB[2]-v[2]*vxB[1],v[2]*vxB[0]-v[0]*vxB[2],v[0]*vxB[1]-v[1]*vxB[0]])# crossproduct
    vxvxB2 = vxvxB2/np.linalg.norm(vxvxB2)
    
    #print "v2, vxB2, vxvxB2 ", v2, vxB2, vxvxB2 
    
    #Backtrafo of efield from shower coord (1,2,3) in xyz (4,5,6) after scaling and/or stretching using the target angles        
    txt1.T[1] = EshowerA.T[0]* v2[0] +EshowerA.T[1]* vxB2[0] + EshowerA.T[2]*vxvxB2[0]
    txt1.T[2] = EshowerA.T[0]* v2[1] +EshowerA.T[1]* vxB2[1] + EshowerA.T[2]*vxvxB2[1]
    txt1.T[3] = EshowerA.T[0]* v2[2] +EshowerA.T[1]* vxB2[2] + EshowerA.T[2]*vxvxB2[2]

#####


        
        


###################################################
#### stretching of positions
###############################

### Note: if kStrecht !=1 then one has to perform an interpolation of antenna positions to be able to compare the signals target to scaled

    # default parametes of star shape simulation
    #rings = 15
    angles= 8
    rings = len(positions[:,1])/angles
    beta= (360./angles)/180.*np.pi


################################
 # Calculating the new stretched antenna positions in the star shape
    offinz= np.mean(positions[:,2])
    offiny= np.mean(positions[:,1])
    offinx= np.mean(positions[:,0])
    pos= np.zeros([len(positions[:,1]),3])
    
    # rotate into shower coordinates for preparation to get the strechted antenna position to compare to - maybe there is a nicer way...
### original    
    GetUVW = UVWGetter(offinx,offiny,offinz, zen1, az1, phigeo, thetageo)
    for i in np.arange(0,len(positions[:,1])):
        pos[i,:] = GetUVW(positions[i,:], )
    
    ##### oldschool
    #for i in np.arange(0,len(positions[:,1])):
        #pos[i,:] = GetUVW(positions[i,:], offinx,offiny,offinz, zen1, az1, phigeo, thetageo)
        
    ###################################   substitue pos by pos and add if condition here
    if kStretch!=1.:
        r= np.linalg.norm(pos[6]- pos[5]) # first rays are formed
        step=1
        #print r
        ## get the distance between the antenna
        if np.linalg.norm(pos[6]-pos[5]) != 0.5*np.linalg.norm(pos[7]-pos[5]): # first rings are formed
            #print ' not equal', np.linalg.norm(pos[6]-pos[5]), 0.5*np.linalg.norm(pos[7]-pos[5])
            #print "ATTENTION: check whether first rings formed or rays if streching is needed"
            r= np.linalg.norm(pos[16]- pos[8])
            step=8
            #print r

        r=r*kStretch # just STRETCHING the radius of the rings

        # default parametes of star shape simulation
        rings = 15
        angles= 8
        beta= (360./angles)/180.*np.pi

        ## get the stretched positions in teh star star shape pattern
        #pos = np.zeros([rings*angles,3])

        # first rays are formed
        if np.linalg.norm(pos[6]-pos[5]) == 0.5*np.linalg.norm(pos[7]-pos[5]): # first rings are formed
            for n in range(0, angles): # second rings
                for m in range(0, rings): #first rays
                    #if n==0 and m==0:
                        #print " first rays than rings"

                    pos[n*rings+m,1 ]= (m+1)*r* np.cos(n*beta)
                    pos[n*rings+m,2 ]= (m+1)*r* np.sin(n*beta)
                    pos[n*rings+m,0 ]= 0.
                    #if l*angles+ n <8:
                    #print l, n*rings+l,  pos[n*rings+l], pos[n*rings+l]

        if np.linalg.norm(pos[6]-pos[5]) != 0.5*np.linalg.norm(pos[7]-pos[5]): # first rings are formed
                for m in range(0, rings): #sencond rays
                    for n in range(0, angles): # first rings
                        #if n==0 and m==0:
                            #print " first rings than rays"

                        pos[m*angles+n,1 ]= (m+1)*r* np.cos(n*beta) # vxB
                        pos[m*angles+n,2 ]= (m+1)*r* np.sin(n*beta) # vxvxB
                        pos[m*angles+n,0 ]= 0. # along v
    
################# CALCULATION OF NEW POSITION VECTOR

    ### the new/target position vector of the star shape plane
    ### pos_prime= p_prime + v_prime *(d_prime + toX)

    #print('distance to max to plane: ',dist1) ## ==toX: fixed distance in m from Xmax to the plane along the shower axis
    #print(v2) # == v_prime: vector of where target shower goes to
    Xmax_primary = _getXmax(primary, E2, zen2)# approximation based on values from plots for gamma (=e) and protons (=pi) # g/cm2
    #print("xmax value " , Xmax_primary)
    Xmax_height, Xmax_distance = _dist_decay_Xmax(zen2, injh2, Xmax_primary)# 8000.# d_prime: distance from decay point to Xmax

    #print("distance to decay to Xmax in m ", Xmax_distance)
    decay=np.array([0.,0.,injh2]) # ==: decay position as defined in zhaires sim, from DANTOn files

    # new position vector:
    #x2= decay - v2 * (Xmax_distance+ dist1) # to account for going from Zhaires to GRAND conv
    x2= decay + v2 * (Xmax_distance+ dist1) # to account for going from Zhaires to GRAND conv

    
    #print "new array center ", x2
    #print "v ", v2
    #print "decay ", decay, " Xmax_distance ", Xmax_distance , "dist ", dist1, " v2 ", v2, " new_array ", x2 

 ##############################   Backtrafo to XYZ
    ### Now the new 'stretched' positions are calculated in the xyz components, backrotation
    stretch2 = np.zeros([len(pos[:,1]),3])

##### original
    #### backtrafo of positions
    GetXYZ = XYZGetter(x2[0],x2[1],x2[2],zen2, az2,phigeo,thetageo)
    for m in range(0, len(pos[:,1])):
        stretch2[m,:] = GetXYZ(pos[m])
    #stretch2[l]  # the new wanted antenna position after stretching
    

##### old school
        
    #stretch2.T[0] = pos.T[0]* v2[0] +pos.T[1]* vxB2[0] + pos.T[2]*vxvxB2[0] + x2[0]
    #stretch2.T[1] = pos.T[0]* v2[1] +pos.T[1]* vxB2[1] + pos.T[2]*vxvxB2[1] + x2[1]
    #stretch2.T[2] = pos.T[0]* v2[2] +pos.T[1]* vxB2[2] + pos.T[2]*vxvxB2[2] + x2[2]

    

    
    #if l==0:
        #print 'position ref to Xmax: ', x2,' position decay: ', decay, ' shower direction: ', v2, ' distance to Xmax: ', dist1, ' distance between Xmax and plane: ', Xmax_distance
        #print len(x2), l, len(stretch2),


    #print ' scaling done , positions ', stretch2[l]
    return txt1, stretch2[l]



################################################################################

def _scale_run(sim_dir, run, primary, E1, zen1, az1, injh1, dist1,
               E2, zen2, az2, injh2, altitude):
    """Scale the simulated traces of a run to the shower parameters
    """

    # TODO: could the magnetic field be given as an argument?
    phigeo =0*np.pi/180.  # 182.66#; (ie pointing 2.66 degrees East from full
                          # North) # phigeo= 0 from simulations inputfile %
                          # In both EVA & Zhaires, North = magnetic North
    thetageo =(180.-27.05)*np.pi/180. # 152.95*np.pi/180. #27.05*np.pi/180.
                                      # (pointing down)-62.95

#    print ""
#    print "###############################  new plane gets scaled ... ", run

    path = os.path.join(sim_dir, run) # Path to the simulation run which
                                      # shall later be rescaled or so

    def print_parameters(E, dist, zen, az, injh):
        """Formated print of the shower parameter values
        """
        parameters = (
            "Energy: {:.2E} EeV".format(E),
            "distance from Xmax: {:.2f} m".format(dist),
            "zenith: {:.2f} deg".format(zen),
            "azimuth: {:.2f} deg".format(az),
            "height: {:.2f} m".format(injh))
        print ", ".join(parameters)

    # Print the reference parameter values
#    print "# Reference shower", run
#    print_parameters(E1, dist1, zen1, az1, injh1)

    # Print the target parameter values
#    print "Target shower parameters"
#    print_parameters(E2, dist1, zen2, az2, injh2)
#    print "Altitude at decay ", altitude, " m"  # to correct injh2 for it since Earth curvature taken into account

    # Convert the angles from degrees to radians
    zen1, az1, zen2, az2 = map(np.deg2rad, (zen1, az1, zen2, az2))

    ######################################### Stretching the position
    ################### Assume kStretch=1, no change of zen and height
    ##r=r*kStretch # just stretching the radius of the rings

    posfile = os.path.join(path, "antpos.dat")
    positions = np.loadtxt(posfile)
    pos_new = np.full_like(positions, 0.)

    ############################################################################

    # Create the output directory if it doesnt exist
    directory = os.path.join(sim_dir, "scaled_" + run)
    if not os.path.exists(directory):
        os.makedirs(directory)

    end = positions.shape[0]
    #### loop over all antenna positions, outer positions should be skipped
    #### since then the interpolation doesnt work anymore
#    print os.path.join(path, "a0.trace")
    for l in np.arange(0, end):
        # always hand over all need parameters,1 3D pulse, and all antenna
        # positions
        txt3, pos_new[l,:] = _scalingpulse(dist1, E1, az1, zen1, injh1, E2, az2,
                                           zen2, injh2, primary, phigeo,
                                           thetageo, l,  positions, path, altitude)



        ###### Writing to file for later use
        name3 = os.path.join(directory, "a{:}.trace".format(l))
        with open(name3, "w+") as FILE:
            # 1,2,3 Ev,vxb,vxvxb then, 4,5,6 Exyz: since in original raw traces
            # Exyz are also in columns 456
            for i in range(0, len(txt3.T[0])):
                args = (txt3.T[0,i], txt3.T[1,i], txt3.T[2,i], txt3.T[3,i])
                print >> FILE, "%.2f	%1.3e	%1.3e	%1.3e" % args



    #print "scaled traces saved like this: {:}/a0.trace".format(directory)

    # Save as well the posiion file somewhere if you scale the complete star
    # shape pattern
    posfile_new = os.path.join(directory, "antpos.dat")
    with open(posfile_new, "w+") as file_ant:
        for i in xrange(0, end):
            args = (pos_new[i,0], pos_new[i,1], pos_new[i,2])
            print >> file_ant, "%.3f	%.3f	%.3f" % args

    #print end, "antennas scaled, positions saved in:", posfile_new


def scale(sim_dir, primary, energy, zenith, azimuth, injection_height, altitude):
    """Scale all simulated traces to the shower parameters
    """
    # Loop over runs
    steerfile_sim = os.path.join(sim_dir, "MasterIndex")
    with open(steerfile_sim, "r") as f:
        for line in f:
            # Unpack the run settings
            args = line.split()
            run = args[0]
            if not os.path.exists(os.path.join(sim_dir, run)):
                continue
            E1, zen1, az1, injh1, dist1 = map(float, args[2:])

            # Conversion from Aires to GRAND convention
            zen1 = 180. - zen1
            az1 = 180. + az1

            # Scale this run
            _scale_run(sim_dir, run, primary, E1, zen1, az1, injh1, dist1,
                       energy, zenith, azimuth, injection_height, altitude)
