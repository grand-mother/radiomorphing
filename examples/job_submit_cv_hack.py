''' This scripts produces the bash files to submit several jobs reading in json files in a specific folder -> e.g. flat100x100_job1.sh
    In addition a script containing all names of that sripts and submitting them is products (eg. msub flat100x100_job1.sh), called run_RM.sh
    Start the submitting with ./run_RM.sh
    
    -> python job_submit.py
'''

 #MOAB commands            brief explanation
 #--------------------------------------------------------------------------
 #msub                     submits a job and queues it in an input queue.
 #showq                    displays information about active, eligible,
                          #blocked, and/or recently completed jobs.
 #showstart                returns start time of submitted job or resources.
 #checkjob                 displays detailed job state information.
 #mjobctl -c <jobid>                   cancels the job `jobid`.


import os

run="hotspot-150x67km2"
hack=1 # off alpha, beta from json file, if =1: alpha=0, beta=0

### path to where one finds the json files which shall be used
path_json= '/project/fh1-project-huepra/le6232/data/retro/'+run+'/taus/'
#path_json= '/project/fh1-project-huepra/qc8087/radiomorphing/examples/seeds/'+run+'/taus/'


# path to the scripts for radiomorphing
path_script='/project/fh1-project-huepra/qc8087/radiomorphing/examples/'
# python ${path_script}/example_ForHLR.py ${path_json}/events.29939363.1.json 


#get an array with all json files which are in the folder
json_file = [f for f in os.listdir(path_json) if f.endswith('.json')]
#print json_file

## craeting a bash file containing all the submission commands
nname='/project/fh1-project-huepra/qc8087/run_CV_hack.sh'
nfile= open(nname, 'w')
nfile.write('#!/bin/bash\n')

# producing bash files to run on multiple cores
for i in range(0, len(json_file )):

        fname="/project/fh1-project-huepra/qc8087/"+run+'_cv_hack_job'+str(i)+'.sh' # filename of each job
        
        command='msub '+str(fname)+'\n'
        nfile.write(command)
        
        
        file= open(fname, 'w')
        file.write('#!/bin/bash\n')
        file.write('#\n')
        file.write('#---------------------\n')
        file.write('# Scheduler options\n')
        file.write('#---------------------\n')
        file.write('#MSUB -q singlenode\n')
        file.write('#MSUB -l nodes=1:ppn=1\n')
        file.write('#MSUB -l walltime=72:00:00\n')
        file.write('#MSUB -l pmem=1gb\n')
        file.write('#MSUB -m bea\n')
        file.write('#MSUB -M zilles@iap.fr\n')
        file.write('#--------------------- \n')       
        file.write('\n')
        file.write('cd $TMP\n')
        message='echo "Job $MOAB_JOBID runnung in $PWD, started: '  +str(path_script)+'/run_computeVoltage_hack.py '+str(path_json)+'/'+str(json_file[i])+'"\n'
        file.write(message)
        script='python '+str(path_script)+'/run_computeVoltage_hack.py '+str(path_json)+'/'+str(json_file[i]) + ' $TMP '+run + ' '+str(hack)
        file.write(script)
        file.close()
        
        i=i+1
        
nfile.close()