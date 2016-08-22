
"""
The Saturn package
==================

Code to analyse the DREAM ENCODE challenge.
"""

from sklearn.metrics import precision_recall_curve

def recall_at_fdr(y_true, y_score, fdr_cutoff=0.05):
    "From: https://www.synapse.org/#!Synapse:syn6131484/discussion/threadId=549"
    precision, recall, thresholds = precision_recall_curve(y_true, y_score) 
    fdr = 1- precision
    cutoff_index = next(i for i, x in enumerate(fdr) if x <= fdr_cutoff)
    return recall[cutoff_index]
