#!/usr/bin/env python3

'''Analyse DNase data for 1 TF across cell types to predict binding
'''

from __future__ import print_function

import logging
import sys
FORMATSTR = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
logger = logging.getLogger(__name__)
handler = logging.StreamHandler(sys.stderr)
handler.setFormatter(logging.Formatter(FORMATSTR))
logger.addHandler(handler)
logger.setLevel(logging.DEBUG)

import numpy as np
np.random.seed(1337)  # for reproducibility

from keras.preprocessing import sequence
from keras.models import Sequential
from keras.layers import Dense, Dropout, Activation, Lambda
from keras.layers import Embedding
from keras.layers import Convolution1D
from keras.datasets import imdb
from keras import backend as K

from Bio import SeqIO
import pyBigWig
import gzip
import pandas as pd
import pickle
import os

# set parameters:
# max_features = 5000
region_size = 200
batch_size = 32  # was max_len
nb_bases = 4
# embedding_dims = 50
nb_filter = 30
filter_length = 3
hidden_dims = 250
nb_epoch = 2
tf = 'CTCF'
cell = 'H1-hESC'
boundunboundratio = .5
maxtrain = 10000
maxtest = maxtrain

BASECODE = np.array(
    ((1, 0, 0, 0),  # A
     (0, 1, 0, 0),  # G
     (0, 0, 1, 0),  # C
     (0, 0, 0, 1),  # T
     (0, 0, 0, 0)), # N
    dtype=np.float)

def baseasint(b):
    "Encode bases as integers."
    if 'A' == b:
        return 0
    if 'C' == b:
        return 1
    if 'G' == b:
        return 2
    if 'T' == b:
        return 3
    if 'N' == b:
        return 4
    raise(ValueError('Bad base "{}" given'.format(b)))

def encodesequence(seq):
    "Encode a sequence as integers."
    return np.fromiter(map(baseasint, seq.upper()), np.int8, len(seq))

def encodebases(bases):
    """Encode bases as a 4-d vector."""
    return np.take(a=BASECODE, indices=bases, axis=0)

def dataforregion(chrom, start, stop):
    data = np.empty((region_size, nb_bases + 1))
    data[:,0:nb_bases] = encodebases(genome[chrom][start:stop])
    data[:,nb_bases] = dnase.values(chrom, start, stop)
    return data

# Load genome from saved numpy data structure or FASTA otherwise
genomenpz = 'data/hg19.npz'
if os.path.exists(genomenpz):
    logger.info('Loading genome: {}'.format(genomenpz))
    genome = dict(np.load(genomenpz).iteritems())
    # sum(map(len, genome.values()))
else:
    genomefasta = '../Data/annotations/hg19.genome.fa.gz'
    genome = {}
    logger.info('Reading genome FASTA: {}'.format(genomefasta))
    for record in SeqIO.parse(gzip.open(genomefasta, mode='rt'), "fasta"):
        logger.info('Loaded: {}'.format(record.id))
        genome[record.id] = encodesequence(record.seq)
    # sum(map(len, genome.values()))
    # Save numpy genome
    np.savez(genomenpz, **genome)

# Load bigwig
dnasefile = "../Data/DNASE/fold_coverage_wiggles/DNASE.{}.fc.signal.bigwig".format(cell)
logger.info('Loading DNase: {}'.format(dnasefile))
dnase = pyBigWig.open(dnasefile)
dnase.chroms()
dnase.header()
# Mean value in range
dnase.stats("chr1", 1000000, 1003000)
# Max value in range
dnase.stats("chr1", 1000000, 1003000, type="max")
# Max value in chromosome
dnase.stats("chr1", type="max")
# All values in range
dnase.values("chr1", 1000000, 1000030)

# Load ChIP data
chipfile = '../Data/ChIPseq/labels/{}.train.labels.tsv.gz'.format(tf)
logger.info('Loading ChIP-seq: {}'.format(chipfile))
chip = pd.read_table(chipfile)


logger.info('Build model...')
model = Sequential()

# we add a Convolution1D, which will learn nb_filter
# word group filters of size filter_length:
# We have one input for each base (excluding 'N') and one for the DNase information
input_shape = (region_size, nb_bases + 1)
model.add(Convolution1D(input_shape=input_shape,
                        nb_filter=nb_filter,
                        filter_length=filter_length,
                        border_mode='valid',
                        activation='relu',
                        subsample_length=1))

# we use max over time pooling by defining a python function to use
# in a Lambda layer
def max_1d(X):
    return K.max(X, axis=1)

model.add(Lambda(max_1d, output_shape=(nb_filter,)))

# We add a vanilla hidden layer:
model.add(Dense(hidden_dims))
model.add(Dropout(0.2))
model.add(Activation('relu'))

# We project onto a single unit output layer, and squash it with a sigmoid:
model.add(Dense(1))
model.add(Activation('sigmoid'))

logger.info('Compile model...')
model.compile(loss='binary_crossentropy',
              optimizer='adam',
              metrics=['accuracy'])

def balancebound(train, cell):
    "Balance the number of bound/unbound samples in the data."
    # Bound sample indexes
    boundidxs = 'B' == train[cell]
    # How many bound samples?
    num_bound = boundidxs.sum()
    # The limit on the number of unbound samples
    maxunbound = int(num_bound / boundunboundratio)
    # Do we exceed the limit?
    if train.shape[0] - num_bound > maxunbound:
        # Separate into bound and unbound
        bound   = train[   boundidxs.values]
        unbound = train[(~boundidxs).values].sample(maxunbound)
        # Return them mixed together
        result = bound.append(unbound)
        return result.reindex(np.random.permutation(result.index))
    else:
        return train

# Choose train and test data
logger.info('Choose train/test...')
testchrs = ['chr2', 'chr7', 'chr20']
istest = chip['chr'].isin(testchrs)
chiptest = balancebound(chip[istest], cell).sample(maxtest)
chiptrain = balancebound(chip[~istest], cell).sample(maxtrain)
chiptest.shape
chiptrain.shape

def getX(chip):
    result = np.empty((chip.shape[0], region_size, nb_bases + 1))
    for i, (index, row) in enumerate(chip.iterrows()):
        result[i] = dataforregion(row['chr'], row['start'], row['stop'])
    return result

def getY(chip, cell):
    return ('B' == chip[cell]).astype(int)

# Munge data
logger.info('Munge data...')
Xtrain = getX(chiptrain)
Xtest  = getX(chiptest)
Ytrain = getY(chiptrain, cell)
Ytest  = getY(chiptest, cell)

# Fit model
logger.info('Fit model...')
model.fit(Xtrain, Ytrain,
          batch_size=batch_size,
          nb_epoch=nb_epoch,
          validation_data=(Xtest, Ytest))
