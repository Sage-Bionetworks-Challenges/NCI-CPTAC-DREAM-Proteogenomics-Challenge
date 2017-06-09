# Automation of validation and scoring
# Make sure you point to the directory where challenge.py belongs and a log directory must exist for the output
cd ~/NCI-CPTAC-Challenge/docker
#---------------------
#Validate submissions
#---------------------
python docker_challenge.py --acknowledge-receipt --canCancel -u $SYNAPSE_USER -p $SYNAPSE_PASS --send-messages --notifications validate --all >> ~/log/score.log 2>&1

#--------------------
#Score submissions
#--------------------
python docker_challenge.py --canCancel -u $SYNAPSE_USER -p $SYNAPSE_PASS --send-messages --notifications score --all >> ~/log/score.log 2>&1

#--------------------
#Stop submissions
#--------------------
#python docker_challenge.py -u $SYNAPSE_USER -p $SYNAPSE_PASS --send-messages --notifications dockerstop --all >> ~/log/score.log 2>&1
