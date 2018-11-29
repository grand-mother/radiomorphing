
import os
import sys
import numpy as np
from frame import UVWGetter, XYZGetter, get_rotation#, GetUVW
from utils import getCerenkovAngle, load_trace, mag

def _getAirDensity(h):
    ''' Calculation of air density at height h
    '''
    #% Using isothermal Model, as in ZHAireS
    rho_0 = 1.225#; % kg/m3
    M = 0.028966#;  %kg/mol
    g = 9.81#; %ms-2
    T = 288.#; % K
    R = 8.32#;
    rho = rho_0*np.exp(-g*M*h/(R*T))
    return rho

###################################

def _getXmax(primarytype, energy, zen2):
    ''' Calculation of Xmax
    '''
    
    # type of primary (electron or pion, energy in EeV, zen2 (GRAND) in rad
    if primarytype=='electron': # aprroximated by gamma shower
        a=82.5 # g/cm2
        c=342.5 #g/cm2
    if primarytype=='pion': # aprroximated by proton
        a=62.5 # g/cm2
        c=357.5 #g/cm2
    Xmax= a*np.log10(energy*10**6.)+c # energy given as E/EeV* 10**6. to be in TeV
    
 

    return Xmax

def _dist_decay_Xmax(zen2, altitude, Xmax_primary): #zen2: zenith of target shower, altitude counts for xmax position
    ''' Calculation of Xmax height and distance
    '''
    
    #% Using isothermal Model
    rho_0 = 1.225*0.001#; % kg/m3 to 0.001g/cm3: 1g/cm3=1000kg/m3, since X given in g/cm2
    M = 0.028966#;  %kg/mol - 1000g/mol
    g = 9.81#; %ms-2
    T = 288.#; % K
    R = 8.32#; J/K/mol , J=kg m2/s2

    hD=altitude
    Xmax_primary= Xmax_primary# g/cm2 
    gamma=np.pi-zen2 # counterpart of angle
    Re= 6370949 # m, Earth radius
    X=0.
    i=0.
    h=hD
    ai=0
    step=10
    while X< Xmax_primary:
        i=i+1
        ai=i*step #100. #m, distance along the shower axis
        hi= -Re+np.sqrt(Re**2. + ai**2. + hD**2. + 2.*Re*hD - 2*ai*np.cos(gamma) *(Re+hD))
        deltah= abs(h-hi) #(h_i-1 - hi)= delta h
        h=hi # new height, Xmax height in m
        X=X+ rho_0*np.exp(-g*M*hi/(R*T)) * step *100.  # Xmax in g/cm2, density in g/cm3, h: m->100cm, 
    
    return h, ai # Xmax_height in m, Xmax_distance along axis in m

def _scalingfactors(E1, az1, zen1, injh1, E2, az2, zen2, injh2, phigeo,thetageo, altitude, primary):
    ''' calculation of scaling factors 
    '''
    
    #print "altitude scaling ", altitude
    # 2: target shower, 1: generic shower as reference
    #################### Energy scaling
    #% Energy
    kE = E2/E1 # both in 1e18eV

    ############## Azimuth scaling
    #% Azimuth
    Bgeo =  [np.cos(phigeo)*np.sin(thetageo), np.sin(phigeo)*np.sin(thetageo),np.cos(thetageo)]
    vref = [np.cos(az1)*np.sin(zen1), np.sin(az1)*np.sin(zen1), np.cos(zen1)]
    vxB_ref = np.cross(vref,Bgeo)
    vxB_ref = np.linalg.norm(vxB_ref)/(np.linalg.norm(vref)*np.linalg.norm(Bgeo))
    v = [np.cos(az2)*np.sin(zen2), np.sin(az2)*np.sin(zen2), np.cos(zen2)]
    vxB = np.cross(v,Bgeo)
    vxB = np.linalg.norm(vxB)/(np.linalg.norm(v)*np.linalg.norm(Bgeo))
    kAz = vxB/vxB_ref   

    ############### Height+Zenith, distance injection point to xmax
    h_ref=injh1
    h=altitude  # actual altitude wrt sealevel at decay position of target position 
    primary1='electron' #hardcoded: reference primary 
    
    Xmax_primary1 = _getXmax(primary1, E1, zen1)# approximation based on values from plots for gamma (=e) and protons (=pi) # g/cm2
    Xmax_height1, Xmax_distance1 = _dist_decay_Xmax(zen1, injh1, Xmax_primary1)# injh1=altitude for reference shower
    hx_ref = h_ref+Xmax_distance1*np.sin(0.5*np.pi-zen1) #   % Height at reference shower Xmax
    ac_ref = getCerenkovAngle(hx_ref)
    rho_ref = _getAirDensity(hx_ref)

    Xmax_primary2 = _getXmax(primary, E2, zen2)# approximation based on values from plots for gamma (=e) and protons (=pi) # g/cm2
    Xmax_height2, Xmax_distance2 = _dist_decay_Xmax(zen2, altitude, Xmax_primary2) 
    hx = h+Xmax_distance2*np.sin(0.5*np.pi-zen2)#   % Height at target shower Xmax 
    ac = getCerenkovAngle(hx) 
    rho = _getAirDensity(hx) 
    
    kStretch = float(ac)/float(ac_ref)#  % Stretch factor for the antenna grid
    kRho = np.sqrt(rho_ref/rho)
    kHeight = kRho/kStretch
    kAmp=kE*kAz*kHeight


    return kStretch, kE, kAz, kHeight




def _scalingpulse(az1, zen1, az2, zen2,  phigeo, thetageo, l ,path,   kStretch, kE, kAz, kHeight): # hand over parameters from reference shower 1 and target shower 2, the number of the antenna in the star shape you would like to have  and all position of complete starshape (for strechting needed), the path to the folder containing the sim, and for now the frequencies (should be removed if one included the antenna response
    ''' scaling electric fields amplitudes 
    '''

#     #Amplitude scaling factor
#     kAmp=kE*kAz*kHeight

    try:
    ##read in full traces of antenna l: 0:time in ns, 1,2,3:: efield
        txt1 = load_trace(path, l)[:-2,:]
    except IOError:
        print("antenna ID ",str(int(l)), " file doesn't exist")
        sys.exit()

    ############## rotation matrix
    ## Convert efield to shower coordinates to apply the scaling
    #R = get_rotation(zen1, az1, phigeo, thetageo)# original
    #EshowerA = np.dot(txt1[:,1:], R) # original
    
    #### Sciling, kHeight includes 1/kStretch
    #EshowerA.T[0] *= kE * kHeight
    #EshowerA.T[1] *= kE * kAz * kHeight
    #EshowerA.T[2] *= kE * kHeight
    
    ###Backtrafo of efield from shower coord ( after scaling and/or stretching using the target angles
    #Rt = get_rotation(zen2, az2, phigeo, thetageo).T # original
    #v2 = Rt[:,0]# original
    #txt1[:,1:] = np.dot(EshowerA, Rt)# original
    ############## 
    
######## slow way
    inc=thetageo
    az=az1
    zen=zen1
    
    B = np.array([np.cos(phigeo)*np.sin(inc), np.sin(phigeo)*np.sin(inc),np.cos(inc)]) 
    B=B/np.linalg.norm(B)
    v = np.array([np.cos(az)*np.sin(zen),np.sin(az)*np.sin(zen),np.cos(zen)]) 
    v=v/np.linalg.norm(v)
    vxB = np.cross(v,B)
    vxB = vxB/np.linalg.norm(vxB)
    vxvxB = np.cross(v,vxB) 
    vxvxB = vxvxB/np.linalg.norm(vxvxB)
    
    # rotation to showerframe
    EshowerA= np.zeros([len(txt1.T[1]),3])
    EshowerA.T[0]= txt1.T[1]* v[0] +txt1.T[2]*v[1]+ txt1.T[3]*v[2]
    EshowerA.T[1]= txt1.T[1]* vxB[0] +txt1.T[2]*vxB[1]+ txt1.T[3]*vxB[2]
    EshowerA.T[2]= txt1.T[1]* vxvxB[0] +txt1.T[2]*vxvxB[1]+ txt1.T[3]*vxvxB[2]
    
    ### Scaling, kHeight includes 1/kStretch
    EshowerA.T[0] *= kE * kHeight
    EshowerA.T[1] *= kE * kAz * kHeight
    EshowerA.T[2] *= kE * kHeight

    ### define backrotation
#     B = np.array([np.cos(phigeo)*np.sin(inc), np.sin(phigeo)*np.sin(inc),np.cos(inc)]) 
#     B=B/np.linalg.norm(B)
    v2 = np.array([np.cos(az2)*np.sin(zen2),np.sin(az2)*np.sin(zen2),np.cos(zen2)]) 
    v2=v2/np.linalg.norm(v2)
    vxB2 = np.cross(v2,B)
    vxB2 = vxB2/np.linalg.norm(vxB2)
    vxvxB2 = np.cross(v2,vxB2)
    vxvxB2 = vxvxB2/np.linalg.norm(vxvxB2)
    
    
    #Backtrafo of efield from shower coord (1,2,3) after scaling and/or stretching using the target angles, efield now in geographic coordinates    
    txt1.T[1] = EshowerA.T[0]* v2[0] +EshowerA.T[1]* vxB2[0] + EshowerA.T[2]*vxvxB2[0]
    txt1.T[2] = EshowerA.T[0]* v2[1] +EshowerA.T[1]* vxB2[1] + EshowerA.T[2]*vxvxB2[1]
    txt1.T[3] = EshowerA.T[0]* v2[2] +EshowerA.T[1]* vxB2[2] + EshowerA.T[2]*vxvxB2[2]
    
    return txt1



        
        

def _stretchpos(dist1, E1, az1, zen1, injh1, E2, az2, zen2, injh2, primary, phigeo,  thetageo, positions, path, altitude, kStretch): # hand over parameters from reference shower 1 and target shower 2, the number of the antenna in the star shape you would like to have  and all position of complete starshape (for strechting needed), the path to the folder containing the sim, and for now the frequencies (should be removed if one included the antenna response
    ''' stretching of reference positions 
    '''


    # default parametes of star shape simulation
    angles= 8 # 8 rays in start shape pattern - hard coded
    rings = len(positions[:,1])/angles
    beta= (360./angles)/180.*np.pi

#################################

#### NOTE scaled reference position in magnetic coordinates as desired antenna positions
    # Calculating the new stretched antenna positions in the star shape
    offinz= np.mean(positions[:,2])
    offiny= np.mean(positions[:,1])
    offinx= np.mean(positions[:,0])
    pos= np.zeros([len(positions[:,1]),3])

    
    # rotate into shower coordinates for preparation, to get the stretched antenna position to compare to 
    GetUVW = UVWGetter(offinx,offiny,offinz, zen1, az1, phigeo, thetageo)
    #get the radius and scale
    r= np.zeros([len(positions[:,1])])
    for i in np.arange(0,len(positions[:,1])):
        pos[i,:] = GetUVW(positions[i,:], )
        r[i]=mag(pos[i,:])*kStretch
        #print r[i]
      
    # first rays are formed, depends on how star shapepattern was set up
    if np.linalg.norm(pos[6]-pos[5]) == 0.5*np.linalg.norm(pos[7]-pos[5]): 
            for n in range(0, angles): # second rings
                for m in range(0, rings): #first rays
                    pos[n*rings+m,1 ]= r[m*angles]* np.cos(n*beta)
                    pos[n*rings+m,2 ]= r[m*angles]* np.sin(n*beta)
                    pos[n*rings+m,0 ]= 0.
 

    if np.linalg.norm(pos[6]-pos[5]) != 0.5*np.linalg.norm(pos[7]-pos[5]): # first rings are formed
                for m in range(0, rings): #second rays
                    for n in range(0, angles): # first rings

                        pos[m*angles+n,1 ]= r[m*angles]* np.cos(n*beta) # vxB
                        pos[m*angles+n,2 ]= r[m*angles]* np.sin(n*beta) # vxvxB
                        pos[m*angles+n,0 ]= 0. # along v
     
    
################## CALCULATION OF NEW POSITION VECTOR
    ### the new target position vector of the star shape plane
    Xmax_primary = _getXmax(primary, E2, zen2)# approximation based on values from plots for gamma (=e) and protons (=pi) # g/cm2

    Xmax_height, Xmax_distance = _dist_decay_Xmax(zen2, altitude, Xmax_primary)# d_prime: distance from decay point to Xmax, altitude counts here, not injection height
    decay=np.array([0.,0.,injh2]) # ==: decay position as defined in zhaires sim, from DANTON files

    # new position vector:
    v2 = np.array([np.cos(az2)*np.sin(zen2),np.sin(az2)*np.sin(zen2),np.cos(zen2)]) 
    v2=v2/np.linalg.norm(v2)
    x2= decay + v2 * (Xmax_distance+ dist1) 
    


 ##############################   Backtrafo to XYZ
    ### Now the new 'stretched' positions are calculated in the xyz components -> backrotation
    stretch2 = np.zeros([len(pos[:,1]),3])

    #### backtrafo of positions in magnetic coordinates
    GetXYZ = XYZGetter(x2[0],x2[1],x2[2],zen2, az2,phigeo,thetageo)
    for m in range(0, len(pos[:,1])):
        stretch2[m,:] = GetXYZ(pos[m])

    return stretch2



################################################################################

def _scale_run(sim_dir, run, primary, E1, zen1, az1, injh1, dist1,
               E2, zen2, az2, injh2, altitude):
    """Scale the simulated traces of a run to the shower parameters
    """
    # TODO: implement the magnetic field strength as an argument
    ## NOTE: ZHAires and reference shower of Radio morphing are in magnetic coordinates (=> azimuth defined wrt to magnetic North): Antenna positions and azimuth given in simulations are defined by x axis pointing towards magnetic North
    phigeo =2.72*np.pi/180.  # (ie pointing 2.72 degrees East from full North) from simulations inputfile %
    thetageo =(180.-27.05)*np.pi/180. # pointing down

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
    if dist1==4000.: # to reduce printout
        print_parameters(E1, dist1, zen1, az1, injh1)

    # Print the target parameter values
    #    print "Target shower parameters"
        print_parameters(E2, dist1, zen2, az2, injh2)

    # Convert the angles from degrees to radians
    zen1, az1, zen2, az2 = map(np.deg2rad, (zen1, az1, zen2, az2))

    # read-in positions of reference shower
    posfile = os.path.join(path, "antpos.dat")
    positions = np.loadtxt(posfile)
    pos_new = np.zeros(positions.shape)


    ############################################################################

    # Create the output directory dor scaling if it doesnt exist
    directory = os.path.join(sim_dir, "scaled_" + run)
    if not os.path.exists(directory):
        os.makedirs(directory)
    
    ### get the scaling factors
    kStretch, kE, kAz, kHeight = _scalingfactors(E1, az1, zen1, injh1, E2, az2, zen2, injh2, phigeo, thetageo, altitude, primary)#, xmax)
        
    end = positions.shape[0]
    #### loop over all antenna positions
    for l in np.arange(0, end):
        # always hand over all antenna positions
        txt3 = _scalingpulse(az1, zen1, az2, zen2,  phigeo, thetageo, l ,path,   kStretch, kE, kAz, kHeight)
        
        ###### Writing to file for later use
        name3 = os.path.join(directory, "a{:}.trace".format(l))
        with open(name3, "w+") as FILE:
            for i in range(0, len(txt3.T[0])):
                args = (txt3.T[0,i], txt3.T[1,i], txt3.T[2,i], txt3.T[3,i])
                print >> FILE, "%.2f	%1.3e	%1.3e	%1.3e" % args


    # stretching of simulated antenna positions 
    pos_new=_stretchpos(dist1, E1, az1, zen1, injh1, E2, az2, zen2, injh2, primary, phigeo,  thetageo,  positions, path, altitude, kStretch)

    # Save the stretched antenna positions
    posfile_new = os.path.join(directory, "antpos.dat")
    with open(posfile_new, "w+") as file_ant:
        for i in xrange(0, end):
            args = (pos_new[i,0], pos_new[i,1], pos_new[i,2])
            print >> file_ant, "%.3f	%.3f	%.3f" % args



def scale(sim_dir, primary, energy, zenith, azimuth, injection_height, altitude):
    """Scale all simulated traces to the shower parameters
    """
    # Loop over the single star-shape planes in the reference shower
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

