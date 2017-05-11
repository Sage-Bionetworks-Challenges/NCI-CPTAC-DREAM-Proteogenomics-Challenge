import docker
import json
import subprocess
import requests
from synapseclient import File, Folder
import synapseutils as synu
import zipfile
import os
import base64

#Synapse Id of Challenge
CHALLENGE_SYN_ID = "syn8228304"
#Synapse Id of directory that you want the log files to go into
CHALLENGE_LOG_FOLDER = "syn9771357"
CHALLENGE_PREDICTION_FOLDER = "syn8729051"
#These are the volumes that you want to mount onto your docker container
OUTPUT_DIR = '/home/ubuntu/output'
TESTDATA_DIR = '/home/ubuntu/evaluation_data'
TRAINING_DIR = '/home/ubuntu/training_data'
#These are the locations on the docker that you want your mounted volumes to be + permissions in docker (ro, rw)
#It has to be in this format '/output:rw'
MOUNTED_VOLUMES = {OUTPUT_DIR:'/output:rw',
                   TESTDATA_DIR:'/evaluation_data:ro',
                   TRAINING_DIR:'/training_data:ro'}
#All mounted volumes here in a list
ALL_VOLUMES = [OUTPUT_DIR,TESTDATA_DIR,TRAINING_DIR]

## Name of your challenge, defaults to the name of the challenge's project
CHALLENGE_NAME = "NCI-CPTAC DREAM Proteogenomics Challenge"

## Synapse user IDs of the challenge admins who will be notified by email
## about errors in the scoring script
ADMIN_USER_IDS = ['3324230']


config_evaluations = [
#Proteogenomics Subchallenge 1  (8720143)
#Proteogenomics Subchallenge 2a  (8720145)
#Proteogenomics Subchallenge 2b  (8720147)
#Proteogenomics Subchallenge 3a (8720149)
#Proteogenomics Subchallenge 3b  (8720151)

    {
        'id':8720143,
        'score_sh':'/score_sc1.sh'
    },
    {
        'id':8720145,
        'score_sh':'/score_sc2.sh'
    },
    {
        'id':8720147,
        'score_sh':'/score_sc2b.sh'
    },
    {
        'id':8720149,
        'score_sh':'/score_sc3a.sh'
    },
    {
        'id':8720151,
        'score_sh':'/score_sc3b.sh'
    }

]

config_evaluations_map = {ev['id']:ev for ev in config_evaluations}

def getBearerTokenURL(dockerRequestURL, user, password):
    initialReq = requests.get(dockerRequestURL)
    auth_headers = initialReq.headers['Www-Authenticate'].replace('"','').split(",")
    for head in auth_headers:
        if head.startswith("Bearer realm="):
            bearerRealm = head.split('Bearer realm=')[1]
        elif head.startswith('service='):
            service = head.split('service=')[1]
        elif head.startswith('scope='):
            scope = head.split('scope=')[1]
    return("{0}?service={1}&scope={2}".format(bearerRealm,service,scope))

def getAuthToken(dockerRequestURL, user, password):
    bearerTokenURL = getBearerTokenURL(dockerRequestURL, user, password)
    auth = base64.b64encode(user + ":" + password)
    bearerTokenRequest = requests.get(bearerTokenURL,
        headers={'Authorization': 'Basic %s' % auth})
    return(bearerTokenRequest.json()['token'])

def zipdir(path, ziph):
    # ziph is zipfile handle
    for root, dirs, files in os.walk(path):
        for file in files:
            ziph.write(os.path.join(root, file),os.path.join(root, file).replace(path+"/",""))


def dockerValidate(submission, syn, user, password):
    submissionJson = json.loads(submission['entityBundleJSON'])
    assert submissionJson['entity'].get('repositoryName') is not None, "Must submit a docker container"
    dockerRepo = submissionJson['entity']['repositoryName'].replace("docker.synapse.org/","")
    #assert dockerRepo.startswith("docker.synapse.org")
    assert submission.get('dockerDigest') is not None, "Must submit a docker container with a docker sha digest"
    dockerDigest = submission['dockerDigest']
    index_endpoint = 'https://docker.synapse.org'
    #dockerImage = dockerRepo + "@" + dockerDigest

    #Check if docker is able to be pulled
    dockerRequestURL = '{0}/v2/{1}/manifests/{2}'.format(index_endpoint, dockerRepo, dockerDigest)
    token = getAuthToken(dockerRequestURL, user, password)

    resp = requests.get(dockerRequestURL,
                        headers={'Authorization': 'Bearer %s' % token})
    assert resp.status_code == 200, "Docker image + sha digest must exist"
    
    #Must check docker image size
    #Synapse docker registry
    dockerSize = sum([layer['size'] for layer in resp.json()['layers']])
    assert dockerSize/1000000 < 1000, "Docker image must be less than a teribyte"

    #Send email to me if harddrive is full 
    #should be stateless, if there needs to be code changes to the docker agent
    checkPredExist = syn.query('select id from folder where parentId == "%s" and name == "%s"' % (CHALLENGE_PREDICTION_FOLDER, submission.id))
    checkLogExist = syn.query('select id from folder where parentId == "%s" and name == "%s"' % (CHALLENGE_LOG_FOLDER, submission.id))

    if checkPredExist['totalNumberOfResults'] == 0:
        predFolder = syn.store(Folder(submission.id, parent = CHALLENGE_PREDICTION_FOLDER))
        predFolder = predFolder.id
    else:
        predFolder = checkPredExist['results'][0]['folder.id']
    if checkLogExist['totalNumberOfResults'] == 0:
        logFolder = syn.store(Folder(submission.id, parent = CHALLENGE_LOG_FOLDER))
        logFolder = logFolder.id
    else:
        logFolder = checkLogExist['results'][0]['folder.id']      
    for participant in submission.contributors:
        if participant['principalId'] in ADMIN_USER_IDS: 
            access = ['CREATE', 'READ', 'UPDATE', 'DELETE', 'CHANGE_PERMISSIONS', 'MODERATE', 'CHANGE_SETTINGS']
        else:
            access = ['READ']
        #Comment set permissions out if you don't want to allow participants to see the pred files
        #syn.setPermissions(predFolder, principalId = participant['principalId'], accessType = access)
        syn.setPermissions(logFolder, principalId = participant['principalId'], accessType = access)
    return(True, "Your submission has been validated!  As your submission is being scored, please go here: https://www.synapse.org/#!Synapse:%s to check on your log files and resulting prediction files." % predFolder)


def dockerRun(submission, scoring_sh, syn, client):
    allLogs = synu.walk(syn, CHALLENGE_LOG_FOLDER)
    logFolder = allLogs.next()
    logFolderId = [synId for name, synId in logFolder[1] if name == submission.id][0]
    
    allPreds = synu.walk(syn, CHALLENGE_PREDICTION_FOLDER)
    predFolder = allPreds.next()
    predFolderId = [synId for name, synId in predFolder[1] if name == submission.id][0]

    dockerDigest = submission.get('dockerDigest')
    submissionJson = json.loads(submission['entityBundleJSON'])
    dockerRepo = submissionJson['entity']['repositoryName']
    dockerImage = dockerRepo + "@" + dockerDigest

    #Mount volumes
    volumes = {}
    for vol in ALL_VOLUMES:
        volumes[vol] = {'bind': MOUNTED_VOLUMES[vol].split(":")[0], 'mode': MOUNTED_VOLUMES[vol].split(":")[1]}

    # Run docker image
    errors = ""
    try:
        container = client.containers.run(dockerImage, scoring_sh, detach=True,volumes = volumes, network_mode="none")
    except docker.errors.APIError as e:
        container = None
        errors = str(e)
    #Create log file
    LogFileName = submission.id + "_log.txt"
    open(LogFileName,'w+').close()
    #While docker is still running (the docker python client doesn't update status)
    if container is not None:
        errors += container.logs()
        while subprocess.Popen(['docker','inspect','-f','{{.State.Running}}',container.name],stdout = subprocess.PIPE).communicate()[0] == "true\n":
            for line in container.logs(stream=True):
                with open(LogFileName,'a') as logFile:
                    logFile.write(line)
                #Only store log file if > 0bytes
                statinfo = os.stat(LogFileName)
                if statinfo.st_size > 0:
                    ent = File(LogFileName, parent = logFolderId)
                    logs = syn.store(ent)
        #Remove container after being done
        container.remove()
    #if log file wasn't created, there is an issue
    if os.stat(LogFileName) == 0:
        with open(LogFileName, "w") as logs:
            logs.write(errors)
        ent = File(LogFileName, parent = logFolderId)
        logs = syn.store(ent)
    
    #client.images.remove(dockerRepo)
    #Zip up predictions and store it into CHALLENGE_PREDICTIONS_FOLDER
    if len(os.listdir(OUTPUT_DIR)) > 0:
        zipf = zipfile.ZipFile(submission.id + '_predictions.zip', 'w', zipfile.ZIP_DEFLATED)
        zipdir(OUTPUT_DIR, zipf)
        zipf.close()

        ent = File(submission.id + '_predictions.zip', parent = predFolderId)
        predictions = syn.store(ent)
        prediction_synId = predictions.id
        os.system("rm -rf %s/*" % OUTPUT_DIR)
    else:
        prediction_synId = None
    return(prediction_synId, logs.id)



def validate_docker(evaluation, submission, syn, client, user, password):
    """
    Find the right validation function and validate the submission.

    :returns: (True, message) if validated, (False, message) if
              validation fails or throws exception
    """
    config = config_evaluations_map[int(evaluation.id)]

    results = dockerValidate(submission, syn, user, password)
    return(results)

def run_docker(evaluation, submission, syn, client):
    """
    Find the right scoring function and score the submission

    :returns: (score, message) where score is a dict of stats and message
              is text for display to user
    """

    config = config_evaluations_map[int(evaluation.id)]
    prediction_synId, log_synId =  dockerRun(submission,config['score_sh'], syn, client)
    if prediction_synId is not None:
        #Comment top line if you don't want to return the synId of prediction file
        #message = "You can find your prediction file here: https://www.synapse.org/#!Synapse:%s" % prediction_synId
        message = "Your prediction file has been stored, but you will not have access to it."
    else:
        message = "No prediction file generated, please check your log files!"
    return (dict(PREDICTION_FILE=prediction_synId, LOG_FILE = log_synId), message)

