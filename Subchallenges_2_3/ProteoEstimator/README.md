# proteo_estimator

- [Overview](#Overview)
- [Installation](#installation)
- [Usage](#Usage)
- [Arguments](#Arguments)
- [Value](#Value)
- [Note](#Note)

## Overview
We present the first data science competition aiming at predicting protein levels from copy number and transcript levels, as well as phosphorylation levels from protein levels. The winning models outperform standard baseline machine learning methods and simply using the transcript levels as proxy for protein levels with respect to prediction performance on new patient samples.
An in depth analysis revealed associations between the commonly predictive genes and essentiality. We provide all the submitted models to the community for re-use and a web application to explore the result of this challenge to support improved large scale proteogenomic characterization of tumor samples and a better understanding of signaling deregulation.

## Installation
```
pip install proteo_estimator
```

Requires Python3
## Usage
```python
import proteo_estimator as pr

# Subchallenge 2: predicting protein levels from copy number and transcript levels
prediction_file = pr.predict_protein_abundances(
        tumor,
        rna,
        cna,
        output_dir,
        logging=True)
```

## Arguments
  
| Parameter                 | Default       |Type       | Description   |	
| :------------------------ |:-------------:|:-------------|:-------------|
| tumor	       |	           |str	          |Tumor type, options are 'breast' and 'ovarian'
| rna	       |	           |str	          |Absolute file path for rna table. Table must be in TSV format of genes x samples
| cna	       |	           |str	          |Absolute file path for cna table. Table must be in TSV format of genes x samples
| output_dir	       |	           |str	          |Absolute file path for output directory. Prediction table and confidence scores will be saved under this directory as prediction.tsv and confidence.tsv
| logging	       |True	           |bool	          |Print progress to stdout

## Return Value
| Output                 |Type       | Description   |	
| :------------------------|:-------------|:-------------|
| prediction_file	      |str	          |Path to tab-separated file of predicted protein levels in the shape of genes x samples. This file will be saved in the directory passed to the parameter "output_dir" as prediction.tsv

## Note
Please ensure that your docker daemon is running in the background.
All file paths must be absolute.
