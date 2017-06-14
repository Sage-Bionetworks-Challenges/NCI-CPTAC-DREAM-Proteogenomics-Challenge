# Use rpy2 if you have R scoring functions
import rpy2.robjects as robjects
import os
filePath = os.path.join(os.path.dirname(os.path.abspath(__file__)),'scoring_functions.R')
robjects.r("source('%s')" % filePath)
corr_by_row = robjects.r('correlation_by_row')
rmse_by_row = robjects.r('RMSE_by_row')
##-----------------------------------------------------------------------------
##
## challenge specific code and configuration
##
##-----------------------------------------------------------------------------


## A Synapse project will hold the assetts for your challenge. Put its
## synapse ID here, for example
## CHALLENGE_SYN_ID = "syn1234567"
CHALLENGE_SYN_ID = "syn8228300"

## Name of your challenge, defaults to the name of the challenge's project
CHALLENGE_NAME = "NCI-CPTAC Proteogenomics Challenge"

## Synapse user IDs of the challenge admins who will be notified by email
## about errors in the scoring script
ADMIN_USER_IDS = [3324230]


## Each question in your challenge should have an evaluation queue through
## which participants can submit their predictions or models. The queues
## should specify the challenge project as their content source. Queues
## can be created like so:
##   evaluation = syn.store(Evaluation(
##     name="My Challenge Q1",
##     description="Predict all the things!",
##     contentSource="syn1234567"))
## ...and found like this:
##   evaluations = list(syn.getEvaluationByContentSource('syn3375314'))
## Configuring them here as a list will save a round-trip to the server
## every time the script starts and you can link the challenge queues to
## the correct scoring/validation functions.  Predictions will be validated and 

def validate_func(submission, goldstandard_path):
    ##Read in submission (submission.filePath)
    ##Validate submission
    ## MUST USE ASSERTION ERRORS!!! 
    ##eg.
    ## assert os.path.basename(submission.filePath) == "prediction.tsv", "Submission file must be named prediction.tsv"
    ## or raise AssertionError()...
    ## Only assertion errors will be returned to participants, all other errors will be returned to the admin
    return(True,"Passed Validation")

def score1(submission, goldstandard_path):
    ##Read in submission (submission.filePath)
    ##Score against goldstandard
    return(score1, score2, score3)

def score2(submission, goldstandard_path):
    ##Read in submission (submission.filePath)
    ##Score against goldstandard
    corr = corr_by_row(submission.filePath, goldstandard_path)
    rmse = rmse_by_row(submission.filePath, goldstandard_path)
    return(corr, rmse)

def score3(submission, goldstandard_path):
    ##Read in submission (submission.filePath)
    ##Score against goldstandard
    return(score1, score2, score3)

evaluation_queues = [

# Proteogenomics Subchallenge 1 (8720143)
# Proteogenomics Subchallenge 2 (8720145)
# Proteogenomics Subchallenge 3 (8720149)
    {
        'id':8720143,
        'scoring_func':score1
        'validation_func':validate_func
        'goldstandard_path':'pros_ova_proteome_sort_common_gene_6577.txt'
    },
    {
        'id':8720145,
        'scoring_func':score2
        'validation_func':validate_func
        'goldstandard_path':'pros_ova_proteome_sort_common_gene_6577.txt'
    },
    {
        'id':8720149,
        'scoring_func':score3
        'validation_func':validate_func
        'goldstandard_path':'pros_ova_proteome_sort_common_gene_6577.txt'
    },
# Proteogenomics Subchallenge 1 Express Lane (9604716)
# Proteogenomics Subchallenge 2 Express Lane (9604717)
# Proteogenomics Subchallenge 3 Express Lane (9604718)
    {
        'id':9604716,
        'scoring_func':score1
        'validation_func':validate_func
        'goldstandard_path':'pros_ova_proteome_sort_common_gene_6577.txt'

    },
    {
        'id':9604717,
        'scoring_func':score2
        'validation_func':validate_func
        'goldstandard_path':'pros_ova_proteome_sort_common_gene_6577.txt'
    },
    {
        'id':9604718,
        'scoring_func':score3
        'validation_func':validate_func
        'goldstandard_path':'pros_ova_proteome_sort_common_gene_6577.txt'

    }
]
evaluation_queue_by_id = {q['id']:q for q in evaluation_queues}


## define the default set of columns that will make up the leaderboard
LEADERBOARD_COLUMNS = [
    dict(name='objectId',      display_name='ID',      columnType='STRING', maximumSize=20),
    dict(name='userId',        display_name='User',    columnType='STRING', maximumSize=20, renderer='userid'),
    dict(name='entityId',      display_name='Entity',  columnType='STRING', maximumSize=20, renderer='synapseid'),
    dict(name='versionNumber', display_name='Version', columnType='INTEGER'),
    dict(name='name',          display_name='Name',    columnType='STRING', maximumSize=240),
    dict(name='team',          display_name='Team',    columnType='STRING', maximumSize=240)]

## Here we're adding columns for the output of our scoring functions, score,
## rmse and auc to the basic leaderboard information. In general, different
## questions would typically have different scoring metrics.
leaderboard_columns = {}
for q in evaluation_queues:
    leaderboard_columns[q['id']] = LEADERBOARD_COLUMNS + [
        dict(name='score',         display_name='Score',   columnType='DOUBLE'),
        dict(name='rmse',          display_name='RMSE',    columnType='DOUBLE'),
        dict(name='auc',           display_name='AUC',     columnType='DOUBLE')]

## map each evaluation queues to the synapse ID of a table object
## where the table holds a leaderboard for that question
leaderboard_tables = {}


def validate_submission(evaluation, submission):
    """
    Find the right validation function and validate the submission.

    :returns: (True, message) if validated, (False, message) if
              validation fails or throws exception
    """
    config = evaluation_queue_by_id[int(evaluation.id)]
    validated, validation_message = config['validation_func'](submission, config['goldstandard_path'])

    return True, validation_message


def score_submission(evaluation, submission):
    """
    Find the right scoring function and score the submission

    :returns: (score, message) where score is a dict of stats and message
              is text for display to user
    """
    config = evaluation_queue_by_id[int(evaluation.id)]
    score = config['scoring_func'](submission, config['goldstandard_path'])
    #Make sure to round results to 3 or 4 digits
    return (dict(score=round(score[0],4), rmse=score[1], auc=score[2]), "You did fine!")


