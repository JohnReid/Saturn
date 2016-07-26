
## Gene expression

- recommended to use TPM units
- annotations only from ENCODE (GENCODE?) annotation files

## Submissions

- Leaderboard 10 submissions per cell/TF combo
- Final submissions over all genome but only ranked on held-out chromosomes

## Challenges

*Across cell-type*

When across cell-type challenge finished, the data will be released for the held-out cell
types and the *within cell-type challenge* will start. Most published methods for TF binding
prediction try this style. These predictions required for ranking in main challenge.


## Scoring

Expecting:

- AUROC > .9
- AUPRC ~ .7 - .9
- FDR 10% will not get most bound regions
- FDR 50% should have most

Rank per TF/cell combo in each statistic
Sum ranks to summarise that TF/cell combo

## Baseline model

- Bins that don't overlap DNase peaks - negatives
- Bins that do overlap are scored using a linear classifier with a log loss function (SGDClassifier in scikit learn)
    * inputs: Known motifs (log2 motif score), region summary statistics (max, 99%, 95%, 75%, 50%, mean)
    * max DNase fold change across each bin


## Deadline

Will be strict with deadline, not like in past.


