# Use rpy2 if you have R scoring functions
import rpy2.robjects as robjects
import os
import zipfile
import pandas as pd
filePath = os.path.join(os.path.dirname(os.path.abspath(__file__)),'scoring_functions.R')
robjects.r("source('%s')" % filePath)
corr_by_row = robjects.r('correlation_by_row')
rmse_by_row = robjects.r('NRMSE_by_row')
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

def validate_func1(prediction_path, goldstandard_path):
    assert os.path.isfile(prediction_path), "Submission file must be named predictions.tsv"
    assert os.stat(prediction_path).st_size > 0, "Prediction file can't be empty"
    pred = pd.read_csv(prediction_path, sep="\t",index_col=0)
    gold = pd.read_csv(goldstandard_path, sep="\t",index_col=0)

    assert all(~pred.index.duplicated()), "There cannot be any duplicated protein ids"
    assert all(~pred.columns.duplicated()), "There cannot be any duplicated sample ids"
    assert all(pred.index.isin(gold.index)), "All protein ids in the prediction file must be present in the goldstandard file, you have these: %s" % ",".join(set(pred.index[~pred.index.isin(gold.index)].map(str)))
    assert all(gold.columns.isin(pred.columns)), "All sample ids must be predicted for, you are missing: %s" % ",".join(gold.columns[~gold.columns.isin(pred.columns)])
    assert all(~pred.isnull()), "There can't be any null values"

    return(True,"Passed Validation")

def validate_func2(prediction_path, goldstandard_path):
    assert os.path.isfile(prediction_path), "Submission file must be named predictions.tsv"
    assert os.stat(prediction_path).st_size > 0, "Prediction file can't be empty"
    pred = pd.read_csv(prediction_path, sep="\t",index_col=0)
    gold = pd.read_csv(goldstandard_path, sep="\t",index_col=0)
    assert all(~pred.index.duplicated()), "There cannot be any duplicated protein ids"
    assert all(~pred.columns.duplicated()), "There cannot be any duplicated sample ids"
    assert all(pred.index.isin(gold.index)), "All protein ids in the prediction file must be present in the goldstandard file, you have these: %s" % ",".join(set(pred.index[~pred.index.isin(gold.index)].map(str)))
    assert all(gold.columns.isin(pred.columns)), "All sample ids must be predicted for, you are missing: %s" % ",".join(gold.columns[~gold.columns.isin(pred.columns)])
    assert all(~pred.isnull()), "There can't be any null values"

    return(True,"Passed Validation")

def score1(prediction_path, goldstandard_path):
    ##Read in submission (submission.filePath)
    ##Score against goldstandard
    return(score1, score2, score3)

def score2(prediction_path, goldstandard_path):
    ##Read in submission (submission.filePath)
    ##Score against goldstandard
    corr = corr_by_row(prediction_path, goldstandard_path)[0]
    rmse = rmse_by_row(prediction_path, goldstandard_path)[0]
    return(corr, rmse)

def score3(prediction_path, goldstandard_path):
    ##Read in submission (submission.filePath)
    ##Score against goldstandard
    return(score1, score2, score3)

evaluation_queues = [

# Proteogenomics Subchallenge 1 (8720143)
# Proteogenomics Subchallenge 2 (8720145)
# Proteogenomics Subchallenge 3 (8720149)
    {
        'id':8720143,
        'scoring_func':score1,
        'validation_func':validate_func1,
        'goldstandard_path':os.path.join(os.path.dirname(os.path.abspath(__file__)),'data_test_true.txt')
    },
    {
        'id':8720145,
        'scoring_func':score2,
        'validation_func':validate_func2,
        'goldstandard_path':os.path.join(os.path.dirname(os.path.abspath(__file__)),'rescaled_prospective_ova_proteome_filtered_5820.txt')
    },
    {
        'id':8720149,
        'scoring_func':score3,
        'validation_func':validate_func2,
        'goldstandard_path':os.path.join(os.path.dirname(os.path.abspath(__file__)),'Noneyet')
    },
# Proteogenomics Subchallenge 1 Express Lane (9604716)
# Proteogenomics Subchallenge 2 Express Lane (9604717)
# Proteogenomics Subchallenge 3 Express Lane (9604718)
    {
        'id':9604716,
        'scoring_func':None,
        'validation_func':validate_func1,
        'goldstandard_path':os.path.join(os.path.dirname(os.path.abspath(__file__)),'data_test_true.txt')

    },
    {
        'id':9604717,
        'scoring_func':None,
        'validation_func':validate_func2,
        'goldstandard_path':os.path.join(os.path.dirname(os.path.abspath(__file__)),'rescaled_prospective_ova_proteome_filtered_5820.txt')
    },
    {
        'id':9604718,
        'scoring_func':None,
        'validation_func':validate_func2,
        'goldstandard_path':os.path.join(os.path.dirname(os.path.abspath(__file__)),'Noneyet')

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


def validate_submission(syn, evaluation, submission):
    """
    Find the right validation function and validate the submission.

    :returns: (True, message) if validated, (False, message) if
              validation fails or throws exception
    """
    config = evaluation_queue_by_id[int(evaluation.id)]
    validation_func = config['validation_func']

    status = syn.getSubmissionStatus(submission)
    #Filter through to fiind PREDICTION_FILE
    synId = filter(lambda input: input.get('key', None) == "PREDICTION_FILE", status.annotations['stringAnnos'])[0]
    entity = syn.get(synId['value'])

    dirname = os.path.dirname(entity.path)
    zfile = zipfile.ZipFile(entity.path)

    for name in zfile.namelist():
      zfile.extract(name, dirname)

    prediction_path = os.path.join(dirname,'predictions.tsv')
    results, validation_message = validation_func(prediction_path, config['goldstandard_path'])

    return True, validation_message


def score_submission(syn, evaluation, submission):
    """
    Find the right scoring function and score the submission

    :returns: (score, message) where score is a dict of stats and message
              is text for display to user
    """
    config = evaluation_queue_by_id[int(evaluation.id)]
    scoring_func = config['scoring_func']

    status = syn.getSubmissionStatus(submission)
    #Filter through to fiind PREDICTION_FILE
    synId = filter(lambda input: input.get('key', None) == "PREDICTION_FILE", status.annotations['stringAnnos'])[0]
    entity = syn.get(synId['value'])

    dirname = entity.cacheDir
    zfile = zipfile.ZipFile(entity.path)

    for name in zfile.namelist():
      zfile.extract(name, dirname)

    prediction_path = os.path.join(dirname,'predictions.tsv')
    if scoring_func is not None:
        corr, rmse = scoring_func(prediction_path,config['goldstandard_path'])
    #Make sure to round results to 3 or 4 digits
        return(dict(corr=round(corr,4), rmse=round(rmse,4)), "You submission was scored.\ncorr: %s\nrmse: %s" %(round(corr,4),round(rmse,4)))
    else:
        return(dict(), "Your prediction file is in the correct format and can be scored.  Please feel free to submit to the real challenge queues.")



