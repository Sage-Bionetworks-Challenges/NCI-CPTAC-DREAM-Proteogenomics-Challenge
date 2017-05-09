# Automation of validation and scoring
# Make sure you point to the directory where challenge.py belongs and a log directory must exist for the output
cd ./
#---------------------
#Validate submissions
#---------------------
python docker_challenge.py -u "synapse user here" -p "password" --send-messages --notifications validate --all >> log/score.log 2>&1

#--------------------
#Score submissions
#--------------------
python docker_challenge.py -u "synpase user here" -p "password" --send-messages --notifications score --all >> log/score.log 2>&1