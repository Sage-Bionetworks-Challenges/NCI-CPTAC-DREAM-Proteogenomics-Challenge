#!/bin/bash
# Automation of validation and scoring
# Make sure you point to the directory where challenge.py belongs and a log directory must exist for the output
# Set environment variables for synapse username and password
cd ~/NCI-CPTAC-Challenge/docker
#---------------------
#Validate submissions
#---------------------
#Proteogenomics Subchallenge 1 Express (9604716)
python docker_challenge.py --acknowledge-receipt --canCancel -u $SYNAPSE_USER -p $SYNAPSE_PASS --send-messages --notifications validate 9604716 >> ~/log/score.log 2>&1

#--------------------
#Score submissions
#--------------------
python docker_challenge.py --canCancel -u $SYNAPSE_USER -p $SYNAPSE_PASS --send-messages --notifications score 9604716 >> ~/log/score.log 2>&1

#--------------------
#Stop submissions
#--------------------
#python docker_challenge.py -u $SYNAPSE_USER -p $SYNAPSE_PASS --send-messages --notifications dockerstop 9604716 >> ~/log/score.log 2>&1
