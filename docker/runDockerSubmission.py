import synapseclient
from synapseclient.exceptions import *
from synapseclient import Project, Folder, File
import synapseutils as synu

from datetime import datetime, timedelta
from StringIO import StringIO

import argparse
import lock
import json
import os
import sys
import time
import traceback
import docker 
import messages
import subprocess
import zipfile

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
LOG_DIR = os.path.join(SCRIPT_DIR,"log")
if not os.path.exists(LOG_DIR):
    os.mkdir(LOG_DIR)

ADMIN_USER_IDS = ["3324230"]

CHALLENGE_NAME = "MM Challenge"

def zipdir(path, ziph):
    # ziph is zipfile handle
    for root, dirs, files in os.walk(path):
        for file in files:
            ziph.write(os.path.join(root, file),os.path.join(root, file).replace(path+"/",""))

def get_user_name(profile):
    names = []
    if 'firstName' in profile and profile['firstName'] and profile['firstName'].strip():
        names.append(profile['firstName'])
    if 'lastName' in profile and profile['lastName'] and profile['lastName'].strip():
        names.append(profile['lastName'])
    if len(names)==0:
        names.append(profile['userName'])
    return " ".join(names)

def update_single_submission_status(status, add_annotations, force=False):
    """
    This will update a single submission's status
    :param:    Submission status: syn.getSubmissionStatus()

    :param:    Annotations that you want to add in dict or submission status annotations format.
               If dict, all submissions will be added as private submissions
    """
    existingAnnotations = status.get("annotations", dict())
    privateAnnotations = {each['key']:each['value'] for annots in existingAnnotations for each in existingAnnotations[annots] if annots not in ['scopeId','objectId'] and each['isPrivate'] == True}
    publicAnnotations = {each['key']:each['value'] for annots in existingAnnotations for each in existingAnnotations[annots] if annots not in ['scopeId','objectId'] and each['isPrivate'] == False}

    if not synapseclient.annotations.is_submission_status_annotations(add_annotations):
        privateAddedAnnotations = add_annotations
        publicAddedAnnotations = dict()
    else:
        privateAddedAnnotations = {each['key']:each['value'] for annots in add_annotations for each in add_annotations[annots] if annots not in ['scopeId','objectId'] and each['isPrivate'] == True}
        publicAddedAnnotations = {each['key']:each['value'] for annots in add_annotations for each in add_annotations[annots] if annots not in ['scopeId','objectId'] and each['isPrivate'] == False} 
    #If you add a private annotation that appears in the public annotation, it switches 
    if sum([key in publicAddedAnnotations for key in privateAnnotations]) == 0:
        pass
    elif sum([key in publicAddedAnnotations for key in privateAnnotations]) >0 and force:
        privateAnnotations = {key:privateAnnotations[key] for key in privateAnnotations if key not in publicAddedAnnotations}
    else:
        raise ValueError("You are trying to add public annotations that are already part of the existing private annotations: %s.  Either change the annotation key or specify force=True" % ", ".join([key for key in privateAnnotations if key in publicAddedAnnotations]))
    if sum([key in privateAddedAnnotations for key in publicAnnotations]) == 0:
        pass
    elif sum([key in privateAddedAnnotations for key in publicAnnotations])>0 and force:
        publicAnnotations= {key:publicAnnotations[key] for key in publicAnnotations if key not in privateAddedAnnotations}
    else:
        raise ValueError("You are trying to add private annotations that are already part of the existing public annotations: %s.  Either change the annotation key or specify force=True" % ", ".join([key for key in publicAnnotations if key in privateAddedAnnotations]))

    privateAnnotations.update(privateAddedAnnotations)
    publicAnnotations.update(publicAddedAnnotations)

    priv = synapseclient.annotations.to_submission_status_annotations(privateAnnotations, is_private=True)
    pub = synapseclient.annotations.to_submission_status_annotations(publicAnnotations, is_private=False)

    for annotType in ['stringAnnos', 'longAnnos', 'doubleAnnos']:
        if priv.get(annotType) is not None and pub.get(annotType) is not None:
            if pub.get(annotType) is not None:
                priv[annotType].extend(pub[annotType])
            else:
                priv[annotType] = pub[annotType]
        elif priv.get(annotType) is None and pub.get(annotType) is not None:
            priv[annotType] = pub[annotType]

    status.annotations = priv
    return(status)

def findFolder(syn, synapseId, folderName):
    walked = synu.walk(syn, synapseId)
    firstWalked = walked.next()
    folderId = [synId for name, synId in firstWalked[1] if name == folderName][0]
    return(folderId)

def attemptStoreLog(syn, entity):
    try:
        logs = syn.store(entity)
        logSynId = logs.id
    except SynapseHTTPError as e:
        logSynId = None
    return(logSynId)

def dockerRun(syn, client, submission, scoring_sh, challenge_prediction_folder, challenge_log_folder, volumes, output_dir, timeQuota=None):
    logFolderId = findFolder(syn, challenge_log_folder, submission.id)
    # allLogs = synu.walk(syn, challenge_log_folder)
    # logFolder = allLogs.next()
    # logFolderId = [synId for name, synId in logFolder[1] if name == submission.id][0]
    predFolderId = findFolder(syn, challenge_prediction_folder, submission.id)

    # allPreds = synu.walk(syn, challenge_prediction_folder)
    # predFolder = allPreds.next()
    # predFolderId = [synId for name, synId in predFolder[1] if name == submission.id][0]

    dockerDigest = submission.get('dockerDigest')
    submissionJson = json.loads(submission['entityBundleJSON'])
    dockerRepo = submissionJson['entity']['repositoryName']
    dockerImage = dockerRepo + "@" + dockerDigest

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    # Run docker image (Can attach container name if necessary)
    errors = None
    start_run = int(time.time()*1000)
    exceedTimeQuota = False
    try:
        container = client.containers.run(dockerImage, scoring_sh, detach=True,volumes = volumes, name=submission.id + "_t" + str(int(time.time())), network_disabled=True, mem_limit="20g")
    except docker.errors.APIError as e:
        container = None
        errors = str(e) + "\n"

    logFileName = submission.id + "_log.txt"
    logSynId = None
    #Create the logfile
    openLog = open(logFileName,'w').close()

    #While docker is still running (the docker python client doesn't update status)
    #Add sleeps
    if container is not None:
        while subprocess.Popen(['docker','inspect','-f','{{.State.Running}}',container.name],stdout = subprocess.PIPE).communicate()[0] == "true\n":
            logFileText = container.logs()
            with open(logFileName,'w') as logFile:
                logFile.write(logFileText)
            statinfo = os.stat(logFileName)
            #Only store log file if > 0bytes
            if statinfo.st_size > 0 and statinfo.st_size/1000.0 <= 50:
                ent = File(logFileName, parent = logFolderId)
                logSynId = attemptStoreLog(syn, ent)
            time.sleep(60)
            runNow = int(time.time()*1000)
            if timeQuota is not None:
                if (runNow - start_run) > timeQuota:
                    #container.stop()
                    subprocess.call(['docker','stop',container.name])
                    exceedTimeQuota = True

        #Must run again to make sure all the logs are captured
        logFileText = container.logs()
        with open(logFileName,'w') as logFile:
            logFile.write(logFileText)
        statinfo = os.stat(logFileName)
        # #Only store log file if > 0bytes
        # if statinfo.st_size > 0 and statinfo.st_size/1000.0 <= 50:
        #     ent = File(logFileName, parent = logFolderId)
        #     logSynId = attemptStoreLog(syn, ent)

        #Remove container and image after being done
        container.remove()
        try:
            client.images.remove(dockerImage)
        except:
            print("Unable to remove image")

    statinfo = os.stat(logFileName)
    if statinfo.st_size == 0:
        with open(logFileName,'w') as logFile:
            if errors is not None:
                logFile.write(errors)
            else:
                logFile.write("No Logs")
        ent = File(logFileName, parent = logFolderId)
        logSynId = attemptStoreLog(syn, ent)
    elif statinfo.st_size /1000.0 > 50:
        oldLogs = subprocess.check_output(["tail",logFileName])
        with open(logFileName,'w') as logFile:
            logFile.write(oldLogs + "\n\nLogs truncated because it exceeded size limit of 50kb")
    outdirs = os.listdir(output_dir)
    if exceedTimeQuota:
        with open(logFileName,'a') as logFile:
            logFile.write("\n Docker Container killed due to exceeding the timeQuota of %s hours." % str(timeQuota / (1000*60*60.0)))
        os.system("rm -rf %s/*" % output_dir)
    ent = File(logFileName, parent = logFolderId)
    logSynId = attemptStoreLog(syn, ent) 


    # if logSynId is None:
    #     logFile = synu.walk(syn, logFolderId)
    #     logFiles = logFile.next()
    #     if len(logFiles[2]) == 0:
    #         ent = File(logFileName, parent = logFolderId)
    #         logs = syn.store(ent)
    #         logSynId = logs.id
    #     else:
    #         logSynId = logFiles[2][0][1]
        
    #Zip up predictions and store it into CHALLENGE_PREDICTIONS_FOLDER

    if len(os.listdir(output_dir)) > 0:
        exceedTimeQuota = False
        zipf = zipfile.ZipFile(submission.id + '_predictions.zip', 'w', zipfile.ZIP_DEFLATED)
        zipdir(output_dir, zipf)
        zipf.close()
        ent = File(submission.id + '_predictions.zip', parent = predFolderId)
        predictions = syn.store(ent)
        prediction_synId = predictions.id
        os.remove(submission.id + '_predictions.zip')
    else:
        prediction_synId = None
    #Remove log file and prediction file
    os.remove(logFileName)
    os.system("rm -rf %s" % output_dir)
    if prediction_synId is not None:
        message = "Your prediction file has been stored, but you will not have access to it."
    else:
        message = "No prediction file generated, please check your log file: https://www.synapse.org/#!Synapse:%s" % logSynId

    return({"PREDICTION_FILE":prediction_synId, "LOG_FILE":logSynId}, message, exceedTimeQuota)

#Get team/individual name
def getTeam(syn, submission):
    if 'teamId' in submission:
        team = syn.restGET('/team/{id}'.format(id=submission.teamId))
        if 'name' in team:
            return(team['name'])
        else:
            return(submission.teamId)
    elif 'userId' in submission:
        profile = syn.getUserProfile(submission.userId)
        return(get_user_name(profile))
    else:
        return('?')

def run(syn, client, submissionId, configFile, challenge_prediction_folder, challenge_log_folder, output_dir, mountedVolumes, canCancel, timeQuota=None, dry_run=False):
    submission = syn.getSubmission(submissionId)
    status = syn.getSubmissionStatus(submissionId)
    evaluation = syn.getEvaluation(submission.evaluationId)
    logFile = open(os.path.join(LOG_DIR,status['id'] + "_log.txt"),'w')
    # if canCancel:
    #     status.canCancel = True
    # status.status = "EVALUATION_IN_PROGRESS"
    # #Store run time and evaluation in progress
    # startTime = {"RUN_START":int(time.time()*1000)}
    # add_annotations = synapseclient.annotations.to_submission_status_annotations(startTime,is_private=False)
    # status = update_single_submission_status(status, add_annotations, force=True)
    # status = syn.store(status)

    status.status = "INVALID"
    #Create dictionary that mounts the volumes
    volumes = {output_dir: {'bind':'/output','mode':'rw'} }
    for mount in mountedVolumes:
        binds = mount.split(":")
        assert len(binds) == 3, "Mounted volumes must be formated- /full/path:/mountpoint:ro"
        volumes[binds[0]] = {'bind':binds[1], 'mode':binds[2]}
    
    with open(configFile, 'r') as config:
        config_evaluations = json.load(config)['config_evaluations']
    score_sh = [ev['score_sh'] for ev in config_evaluations if ev['id'] == int(evaluation.id)]

    #If submission_info is None, then the code passed
    submission_info = None
    try:
        score, message, exceedTimeQuota = dockerRun(syn, client, submission, score_sh, challenge_prediction_folder, challenge_log_folder, volumes, output_dir, timeQuota)

        logFile.write("scored: %s %s %s %s" % (submission.id, submission.name, submission.userId, str(score)))
        logFile.flush()
        score['team'] = getTeam(syn, submission)
        score['RUN_END'] = int(time.time()*1000)
        if exceedTimeQuota:
            score['FAILURE_REASON'] = "Exceeded Time Quota of %s hours" % str(timeQuota / (1000*60*60.0))

        add_annotations = synapseclient.annotations.to_submission_status_annotations(score,is_private=False)
        status = update_single_submission_status(status, add_annotations, force=True)
        if score['PREDICTION_FILE'] is None:
            status.status = "INVALID"
        else:
            status.status = "ACCEPTED"

    except Exception as ex1:
        logFile.write('\n\nError scoring submission %s %s:\n' % (submission.name, submission.id))
        st = StringIO()
        traceback.print_exc(file=st)
        message = st.getvalue()
        logFile.write(message)
        logFile.flush()

        if ADMIN_USER_IDS:
            submission_info = "submission id: %s\nsubmission name: %s\nsubmitted by user id: %s\n\n" % (submission.id, submission.name, submission.userId)
            messages.error_notification(userIds=ADMIN_USER_IDS, message=submission_info+st.getvalue(),queue_name=evaluation.name)

    if not dry_run:
        status = syn.store(status)

    ## send message AFTER storing status to ensure we don't get repeat messages
    profile = syn.getUserProfile(submission.userId)
    if status.status == 'ACCEPTED':
        messages.scoring_succeeded(
            userIds=[submission.userId],
            message=message,
            username=get_user_name(profile),
            queue_name=evaluation.name,
            submission_name=submission.name,
            submission_id=submission.id)
    elif submission_info is None:
        messages.scoring_error(
            userIds=[submission.userId],
            message=message,
            username=get_user_name(profile),
            queue_name=evaluation.name,
            submission_name=submission.name,
            submission_id=submission.id)

def command_run(args):
    run(args.syn, args.client, args.submissionId, args.configFile, args.challengePredFolder, args.challengeLogFolder, args.outputDir, args.mountedVolumes, args.canCancel, timeQuota=args.timeQuota, dry_run=args.dry_run)

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("submissionId", metavar="Submission Id")
    parser.add_argument("--challengePredFolder",required=True)
    parser.add_argument("--challengeLogFolder",required=True)
    parser.add_argument("--outputDir", required=True)
    parser.add_argument("--mountedVolumes", nargs="*", required=True)
    parser.add_argument("--configFile", required=True)
    parser.add_argument("--timeQuota", help="Time quota in milliseconds", type=int)
    #Has default values
    parser.add_argument("-u", "--user", help="UserName", default=None)
    parser.add_argument("-p", "--password", help="Password", default=None)
    parser.add_argument("--notifications", help="Send error notifications to challenge admins", action="store_true", default=False)
    parser.add_argument("--send-messages", help="Send validation and scoring messages to participants", action="store_true", default=False)
    parser.add_argument("--acknowledge-receipt", help="Send confirmation message on passing validation to participants", action="store_true", default=False)
    parser.add_argument("--dry-run", help="Perform the requested command without updating anything in Synapse", action="store_true", default=False)
    parser.add_argument("--debug", help="Show verbose error output from Synapse API calls", action="store_true", default=False)
    parser.add_argument("--canCancel", action="store_true", default=False)

    #Test run
    #python runDockerSubmission.py 9636069 --challengePredFolder syn7998461 --challengeLogFolder syn9974718 --configFile config.json --mountedVolumes /home/ubuntu/sc2/Celgene-Multiple-Myeloma-Challenge/docker_agent/test-data:/test-data:ro /.synapseCache:/.synapseCache:ro -u $SYNAPSE_USER -p $SYNAPSE_PASS --outputDir /home/ubuntu/sc2/Celgene-Multiple-Myeloma-Challenge/docker_agent/9636069 --send-messages --notifications   
    args = parser.parse_args()

    print "\n" * 2, "=" * 75
    print datetime.utcnow().isoformat()

    # Acquire lock, don't run two scoring scripts at once
    try:
        submission_lock = lock.acquire_lock_or_fail(args.submissionId, max_age=timedelta(hours=9000))
    except lock.LockedException:
        print u"Is the scoring script already running? Can't acquire lock."
        # can't acquire lock, so return error code 75 which is a
        # temporary error according to /usr/include/sysexits.h
        return 75

    try:
        syn = synapseclient.Synapse(debug=args.debug)
        if not args.user:
            args.user = os.environ.get('SYNAPSE_USER', None)
        if not args.password:
            args.password = os.environ.get('SYNAPSE_PASSWORD', None)
        syn.login(email=args.user, password=args.password)
        #Add client into arguments
        client = docker.from_env()
        client.login(args.user, args.password, registry="http://docker.synapse.org")
        #Add syn and client into arguments
        args.syn = syn
        args.client = client

        ## initialize messages
        messages.syn = syn
        messages.dry_run = args.dry_run
        messages.send_messages = args.send_messages
        messages.send_notifications = args.notifications
        messages.acknowledge_receipt = args.acknowledge_receipt
        command_run(args)

    except Exception as ex1:
        sys.stderr.write('Error in scoring script:\n')
        st = StringIO()
        traceback.print_exc(file=st)
        sys.stderr.write(st.getvalue())
        sys.stderr.write('\n')
        if ADMIN_USER_IDS:
            messages.error_notification(userIds=ADMIN_USER_IDS, message=st.getvalue(), queue_name=CHALLENGE_NAME)

    finally:
        submission_lock.release()

    print "\ndone: ", datetime.utcnow().isoformat()
    print "=" * 75, "\n" * 2


if __name__ == '__main__':
    main()