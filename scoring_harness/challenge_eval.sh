# Automation of validation and scoring
# Make sure you point to the directory where challenge.py belongs and a log directory must exist for the output
cd ./
#---------------------
#Validate submissions
#---------------------
#Remove --send-messages to do rescoring without sending emails to participants
python challenge.py -u "synapse user here" --send-messages --notifications validate --all >> log/score.log 2>&1

#--------------------
#Score submissions
#--------------------
python challenge.py -u "synpase user here" --send-messages --notifications score --all >> log/score.log 2>&1