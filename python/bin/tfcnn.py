#!/usr/bin/env python3

'''Analyse DNase data for 1 TF across cell types to predict binding
'''

from __future__ import print_function
np.random.seed(1337)  # for reproducibility

from keras.preprocessing import sequence
from keras.models import Sequential
from keras.layers import Dense, Dropout, Activation, Lambda
from keras.layers import Embedding
from keras.layers import Convolution1D
from keras.datasets import imdb
from keras import backend as K

import numpy as np
from Bio import SeqIO
import pyBigWig
import gzip
import pandas as pd
import pickle

# set parameters:
# max_features = 5000
region_size = 200
batch_size = 32  # was max_len
nb_bases = 5  # Include 'N'
# embedding_dims = 50
nb_filter = 30
filter_length = 3
hidden_dims = 250
nb_epoch = 2
# We have 5 inputs for each base and one for the DNase information
input_shape = (region_size, nb_bases + 1)

print('Build model...')
model = Sequential()

# we add a Convolution1D, which will learn nb_filter
# word group filters of size filter_length:
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

model.compile(loss='binary_crossentropy',
              optimizer='adam',
              metrics=['accuracy'])

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

# Load genome from saved numpy data structure or FASTA otherwise
genomenpz = 'data/hg19.npz'
if os.path.exists(genomenpz):
    genome = dict(np.load(genomenpz).iteritems())
    # sum(map(len, genome.values()))
else:
    genomefasta = '../Data/annotations/hg19.genome.fa.gz'
    genome = {}
    for record in SeqIO.parse(gzip.open(genomefasta, mode='rt'), "fasta"):
        print('Loaded: {}'.format(record.id))
        genome[record.id] = encodesequence(record.seq)
    # sum(map(len, genome.values()))
    # Save numpy genome
    np.savez(genomenpz, **genome)

# Load bigwig
dnase = pyBigWig.open("../Data/DNASE/fold_coverage_wiggles/DNASE.H1-hESC.fc.signal.bigwig")
dnase.chroms()
dnase.header()
# Mean value in range
dnase.stats("chr1", 1000000, 1003000)
# Max value in range
dnase.stats("chr1", 1000000, 1003000, type="max")
# Max value in chromosome
dnase.stats("chr1", type="max")
# All values in range
dnase.values("chr1", 1000000, 1003000)
dnase.values("chr1", 1000000, 1000030)
dnase.values("chr1", 0, 3)

# Load ChIP data
chip = pd.read_table('../Data/ChIPseq/labels/CTCF.train.labels.tsv.gz')
state = 'B'
cell = 'H1-hESC'
num_sample = 10
def samplebound(chip, state, cell, num_sample):
    "Sample rows from chip data which match the cell's state"
    idxs = np.argwhere(state == chip[cell])[:,0]
    sampledidxs = np.random.choice(idxs, size=num_sample, replace=False)
    sampledidxs
    return chip.iloc[sampledidxs]
bound = samplebound(chip, state, cell, num_sample)
bound

# Munge data
BASECODE = np.array(
    ((1, 0, 0, 0),  # A
     (0, 1, 0, 0),  # G
     (0, 0, 1, 0),  # C
     (0, 0, 0, 1),  # T
     (0, 0, 0, 0)), # N
    dtype=np.float)
def encodebases(bases):
    """Encode bases as a 4-d vector."""
    return np.take(a=BASECODE, indices=bases, axis=0)
encodebases([0,1,2,2,3,4])
tmp = encodebases(genome['chr20'])
tmp[1000000:1000100]

def dataforregion(chrom, start, stop):
    data = np.empty((200, nb_bases + 1))
    data[:,0:nb_bases] = genome[chrom][start:stop]
    data[:,nb_bases] = dnase.values(chrom, start, stop)
    return data
dataforregion('chr1')

# Fit model
model.fit(X_train, y_train,
          batch_size=batch_size,
          nb_epoch=nb_epoch,
          validation_data=(X_test, y_test))
