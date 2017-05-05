import docker
import json
import subprocess
from synapseclient import File
import zipfile
import os
import base64

#Synapse Id of Challenge
CHALLENGE_SYN_ID = "syn1235"
#Synapse Id of directory that you want the log files to go into
CHALLENGE_LOG_FOLDER = "syn8729051"
CHALLENGE_PREDICTION_FOLDER = "syn8729051"
#These are the volumes that you want to mount onto your docker container
OUTPUT_DIR = '/home/ubuntu/output'
TESTDATA_DIR = '/home/ubuntu/test-data'
#These are the locations on the docker that you want your mounted volumes to be + permissions in docker (ro, rw)
#It has to be in this format '/output:rw'
MOUNTED_VOLUMES = {OUTPUT_DIR:'/output:rw',
				   TESTDATA_DIR:'/test-data:ro'}
#All mounted volumes here in a list
ALL_VOLUMES = [OUTPUT_DIR,TESTDATA_DIR]

## Name of your challenge, defaults to the name of the challenge's project
CHALLENGE_NAME = "Example Synapse Challenge"

## Synapse user IDs of the challenge admins who will be notified by email
## about errors in the scoring script
ADMIN_USER_IDS = ['123234']


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
    auth_headers = initialReq.headers['Www-Authenticate'].split(",")
    auth_headers = [i.split("=")[1].replace('"','') for i in auth_headers]
    required = ['Bearer realm=','service=','scope=']
    for head in auth_headers:
        if head.startswith("Bearer realm="):
            bearerRealm = head.split('Bearer realm=')[1]
        elif head.startswith('service='):
            service = head.split('service=')[1]
        elif head.startswith('scope='):
            scope = head.split('scope=')[1]
    return("{0}?service={1}&scope={2}".format(bearerRealm,service,scope))

def getAuthToken(dockerRequestURL, user, password):
    bearerTokenURL = getBearerTokenURL(dockerRequestURL)
    auth = base64.b64encode(user:password)
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
    dockerRepo = submissionJson['entity']['repositoryName']
    #assert dockerRepo.startswith("docker.synapse.org")
    assert submission.get('dockerDigest') is not None, "Must submit a docker container with a docker sha digest"
    dockerDigest = submission['dockerDigest']
    index_endpoint = 'https://docker.synapse.org'
    #dockerImage = dockerRepo + "@" + dockerDigest

    #Check if docker is able to be pulled
    dockerRequestURL = '{0}/v2/{1}/manifests/{2}'.format(index_endpoint, dockerRepo, dockerDigest)
    token = getAuthToken(dockerRequestURL, user, password)

    resp = requests.get('{0}/v2/{1}/manifests/{2}'.format(index_endpoint, dockerRepo, dockerDigest),
                        headers={'Authorization': 'Bearer %s' % token})

    assert resp.status_code == 200, "Docker image + sha digest must exist"

    #Must check docker image size
    #Send email to me if harddrive is full 
    #should be stateless, if there needs to be code changes to the docker agent

    checkExist = syn.query('select id from folder where parentId == "%s" and name == "%s"' % (CHALLENGE_LOG_PREDICTION_FOLDER, submission.id))
    if checkExist['totalNumberOfResults'] == 0:
        predFolder = syn.store(Folder(submission.id, parent = CHALLENGE_LOG_PREDICTION_FOLDER))
        predFolder = predFolder.id
    else:
        predFolder = checkExist['results'][0]['folder.id']
    for participant in submission.contributors:
        if participant['principalId'] in ADMIN_USER_IDS: 
            access = ['CREATE', 'READ', 'UPDATE', 'DELETE', 'CHANGE_PERMISSIONS', 'MODERATE', 'CHANGE_SETTINGS']
        else:
            access = ['READ']
        syn.setPermissions(predFolder, principalId = participant['principalId'], accessType = access)

    return(True, "Your submission has been validated!  As your submission is being scored, please go here: https://www.synapse.org/#!Synapse:%s to check on your log files and resulting prediction files." % predFolder)


def dockerRun(submission, scoring_sh, syn, client):
    allSubmissions = synu.walk(syn, CHALLENGE_LOG_PREDICTION_FOLDER)
    predFolder = allSubmissions.next()
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
	container = client.containers.run(dockerRepo, 'bash /score_sc1.sh', detach=True,volumes = volumes, network_disabled=True)
	
	#Create log file
	LogFileName = submission.id + "_log.txt"
	open(LogFileName,'w+').close()
	
	#While docker is still running (the docker python client doesn't update status)
	while subprocess.Popen(['docker','inspect','-f','{{.State.Running}}',container.name],stdout = subprocess.PIPE).communicate()[0] == "true\n":
		for line in container.logs(stream=True):
			with open(LogFileName,'a') as logFile:
				logFile.write(line)
			#Only store log file if > 0bytes
			statinfo = os.stat(LogFileName)
			if statinfo.st_size > 0:
				ent = File(LogFileName, parent = CHALLENGE_LOG_FOLDER)
				logs = syn.store(ent)

	#Remove container after being done
    container.remove()
    client.images.remove(dockerImage)
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


def run_docker(evaluation, submission):
	"""
	Find the right scoring function and score the submission

	:returns: (score, message) where score is a dict of stats and message
			  is text for display to user
	"""

    config = config_evaluations_map[int(evaluation.id)]
    prediction_synId, log_synId =  dockerRun(submission,config['score_sh'], syn, client)
    if prediction_synId is not None:
        message = "You can find your prediction file here: https://www.synapse.org/#!Synapse:%s" % prediction_synId
    else:
        message = "No prediction file generated, please check your log files!"
    return (dict(PREDICTION_FILE=prediction_synId, LOG_FILE = log_synId), message)

