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
from keras import backend as K
from keras import layers as la

from Bio import SeqIO
import pyBigWig
import gzip
import pandas as pd
import pickle
import os
import functools
import operator
import scipy.stats as st
from prg import prg

# set parameters:
# max_features = 5000
region_size = 200
batch_size = 32  # was max_len
nb_bases = 4
# embedding_dims = 50
nb_filter = 50
filter_length = 6 * nb_bases
hidden_dims = 250
nb_epoch = 1000
tf = 'CTCF'
cell = 'H1-hESC'
boundunboundratio = .5
maxtrain = 10000
maxtest = maxtrain

def noise_for(x, scale):
    return np.random.normal(scale=.1, size=x.shape)


###################
# Load data section
###################
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


##################################
# Split training/test data section
##################################
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
        # Return them, no need to shuffle as Keras will do this for us
        return bound.append(unbound)
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
    """Encode bases as a 1-d vector with 4 entries for each base."""
    return np.take(a=BASECODE, indices=bases, axis=0).flatten()

def getseqinput(chip):
    result = np.empty((chip.shape[0], region_size * nb_bases, 1))
    for i, (index, row) in enumerate(chip.iterrows()):
        result[i,:,0] = encodebases(genome[row['chr']][row['start']:row['stop']])
    return result

def getdnaseinput(chip):
    result = np.empty((chip.shape[0], region_size, 1))
    for i, (index, row) in enumerate(chip.iterrows()):
        result[i,:,0] = dnase.values(row['chr'], row['start'], row['stop'])
        # if np.isnan(result[i]).sum() != 0:
        #     raise ValueError('Found NaN in DNase data: %s:%d-%d' %
        #             (row['chr'], row['start'], row['stop']))
    return np.nan_to_num(result)

def getchipoutput(chip, cell):
    return np.where(('B' == chip[cell]).values, 1., 0.)

# Munge data
logger.info('Munge data...')
seqinputtrain = getseqinput(chiptrain)
seqinputtest  = getseqinput(chiptest)
dnaseinputtrain = getdnaseinput(chiptrain)
dnaseinputtest  = getdnaseinput(chiptest)
outtrain = getchipoutput(chiptrain, cell)
outtest  = getchipoutput(chiptest , cell)


###############
# Model section
###############
#
# Sequence branch
#
logger.info('Build model...')
seqbranch = Sequential()
#
# we add a Convolution1D, which will learn nb_filter
# word group filters of size filter_length:
# We have one input for each base (excluding 'N') and one for the DNase information
input_shape = (region_size * nb_bases,1)
seqbranch.add(la.Convolution1D(
    input_shape=input_shape,
    nb_filter=nb_filter,
    filter_length=filter_length,
    border_mode='valid',
    #activation='relu',
    activation='tanh',
    subsample_length=1))
seqbranch.output_shape
#
# Average each convolutional feature over the whole region
seqbranch.add(la.AveragePooling1D(pool_length=seqbranch.output_shape[1]))

#
# DNase branch
#
dnasebranch = Sequential()
dnasebranch.add(la.Convolution1D(
    input_shape=(region_size, 1),
    nb_filter=5,
    filter_length=25,
    border_mode='valid',
    #activation='relu',
    activation='tanh',
    subsample_length=16))
dnasebranch.output_shape
#
# Common layer
#
model = dnasebranch
# We add a vanilla hidden layer:
model.add(la.Flatten())
model.add(la.Dense(5))
model.add(la.Dropout(0.2))
model.add(la.Activation('relu'))
#
# We project onto a single unit output layer, and squash it with a sigmoid:
model.add(la.Dense(1))
model.add(la.Activation('sigmoid'))
#
logger.info('Compile model with %d parameters...', model.count_params())
model.compile(
    loss='binary_crossentropy',
    optimizer='adam',
    # optimizer='rmsprop',
    metrics=['accuracy'])
#
#
###########
# Fit model
###########
logger.info('Fit model...')
hist = model.fit(
    dnaseinputtrain, outtrain,
    batch_size=batch_size,
    nb_epoch=5000,
    # nb_epoch=10,
    validation_data=(dnaseinputtest, outtest))
logger.info('Done')
print(hist.history)
hist.history['val_loss']
hist.history['val_acc']

#################
# Analyse results
#################
scores = model.predict(dnaseinputtest)
prg_curve = prg.create_prg_curve(outtest, scores.flatten())
auprg = prg.calc_auprg(prg_curve)
logger.info('AUPRG = %f', auprg)
fig = prg.plot_prg(prg_curve)
fig.savefig('AUPRG.png')
