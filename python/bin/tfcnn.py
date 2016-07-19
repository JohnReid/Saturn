#!/usr/bin/env python3

'''Analyse DNase data for 1 TF across cell types to predict binding
'''

from __future__ import print_function
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

print('Build model...')
model = Sequential()

# We have 5 inputs for each base and one for the DNase information
model.add(Dense(nb_bases + 1, input_shape=(region_size,)))

# we add a Convolution1D, which will learn nb_filter
# word group filters of size filter_length:
model.add(Convolution1D(nb_filter=nb_filter,
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

print('Loading data...')
(X_train, y_train), (X_test, y_test) = imdb.load_data(nb_words=max_features,
                                                      test_split=0.2)
print(len(X_train), 'train sequences')
print(len(X_test), 'test sequences')

print('Pad sequences (samples x time)')
X_train = sequence.pad_sequences(X_train, maxlen=maxlen)
X_test = sequence.pad_sequences(X_test, maxlen=maxlen)
print('X_train shape:', X_train.shape)
print('X_test shape:', X_test.shape)

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

def encodesequence(seq):
    "Encode a sequence as integers."
    return numpy.fromiter(map(baseasint, seq.upper()), np.int8, len(seq))

# Load genome
genomefasta = '../Data/annotations/hg19.genome.fa.gz'
genome = {}
with gzip.open(genomefasta) as handle:
    for record in SeqIO.parse(handle, "fasta"):
        genome[record.id] = encodesequence(record.seq)
genome

# Load bigwig
bw = pyBigWig.open("../Data/DNASE/fold_coverage_wiggles/DNASE.H1-hESC.fc.signal.bigwig")
bw.chroms()
bw.header()
# Mean value in range
bw.stats("chr1", 1000000, 1003000)
# Max value in range
bw.stats("chr1", 1000000, 1003000, type="max")
# Max value in chromosome
bw.stats("chr1", type="max")
# All values in range
bw.values("chr1", 1000000, 1003000)
bw.values("chr1", 1000000, 1000030)
bw.values("chr1", 0, 3)

# Load binding data
import pandas as pd
binding = pd.read_table('../Data/ChIPseq/labels/CTCF.train.labels.tsv.gz')

# Fit model
model.fit(X_train, y_train,
          batch_size=batch_size,
          nb_epoch=nb_epoch,
          validation_data=(X_test, y_test))
