#!/usr/bin/env Rscript
#
# Make features from DNAshape bigwig files
#
#!
#! sbatch directives begin here ###############################
#!
#SBATCH -p mrc-bsu-sand                 # Partition
#SBATCH -J FEATSHAPE                    # Name of the job
#SBATCH -A MRC-BSU-SL2                  # Which project should be charged
#SBATCH --nodes=1                       # How many whole nodes should be allocated?
#SBATCH --ntasks=1                      # How many (MPI) tasks will there be in total? (<= nodes*16)
#SBATCH --mem=15360                     # How many MB each node is allocated
#SBATCH --time=06:00:00                 # How much wallclock time will be required?
#SBATCH -o feat-shape/feat-shape-%j.out # stdout
#SBATCH -e feat-shape/feat-shape-%j.out # stderr
#SBATCH --mail-type=FAIL                # What types of email messages do you wish to receive?
##SBATCH --no-requeue                   # Uncomment this to prevent the job from being requeued (e.g. if
"Usage: feature-shape.R [--wellington] FEATURE BIGWIG..." -> doc


#
# Set warnings as errors
#
options(warn = 2)


#
# Load packages
#
devtools::load_all()
library(Saturn)
library(stringr)
library(rtracklayer)


#
# Parse options
#
# .args <- "--wellington DNASE ../Data/DNASE/fold_coverage_wiggles/DNASE.SK-N-SH.fc.signal.bigwig"
if (! exists(".args")) .args <- commandArgs(TRUE)
opts <- docopt::docopt(doc, args = .args)
print(opts)
feature <- opts[['FEATURE']]
bigwig.file.paths <- opts[['BIGWIG']]
wellington <- opts$wellington


#
# Set up
#
for (big.file.path in bigwig.file.paths) {
  message('BigWig file: ', bigwig.file.path)
  stopifnot(file.exists(bigwig.file.path))
}


#' Cached function to load bigwig
#'
bw.file <- memoise::memoise(function(bigwig.file.path) BigWigFile(bigwig.file.path))


#' Generate feature for each set of ranges
#'
ranges.test.pp <- GNCList(ranges.test())
calc.feature <- function(hits.gr, score.name = 'score') {
  overlaps <- as.data.frame(findOverlaps(ranges.test.pp, hits.gr, ignore.strand = TRUE))
  if (! nrow(overlaps)) {
    Rle(0, length(ranges.test.pp))
  } else {
    overlaps$value <- mcols(hits.gr[overlaps$subjectHits])[[score.name]]
    with(
      aggregate(overlaps %>% select(-subjectHits), list(overlaps$queryHits), max),
      Rle.from.sparse(length(ranges.test.pp), queryHits, value))
  }
}


#' Subset ranges by footprints
#'
subsetByFootprints <- function(gr, cell) {
  footprints <- load.wellington(cell)
  subsetByOverlaps(gr, footprints)
}


#' Save features as sparse Matrix
#'
save.features <- function(features, features.file.name) {
  message('Saving features: ', features.file.name)
  dir.create(dirname(features.file.name), showWarnings = FALSE)
  saveRDS(as(do.call(S4Vectors::DataFrame, features), features.file.name, "Matrix"))
}


#
# Are we subsetting by Wellington footprints?
#
if (wellington) {
  for (cell in cells) {
    message('Calculating features in: ', cell)
    #
    # If a TF is specified name the feature file with it
    #
    well.tag <- str_c(scan.tag, 'Well')
    if (is.null(tf)) {
      features.file.name <- feature.cell.file.name(well.tag, cell)
    } else {
      features.file.name <- feature.tf.cell.file.name(well.tag, tf, cell)
    }
    if (file.exists(features.file.name)) {
      message('Features already exist, not recreating: ', features.file.name)
    } else {
      #
      # Calculate and save features (subsetted by footprints)
      message('Subsetting by footprints: ', cell)
      features <-
        lapply(
          lapply(get.scan(), functional::Curry(subsetByFootprints, cell = cell)),
          calc.feature)
      names(features) <- str_c(names(features), '.Well')
      save.features(features, features.file.name)
    }
  }
}


#
# Feature name from file name
#
feature.from.file <- function(file.name) {
  basename(file.name) %>%
    str_replace(regex('.bw$'    , ignore_case = TRUE), '') %>%
    str_replace(regex('.bigwig$', ignore_case = TRUE), '') %>%
    str_replace(regex('.wig$'   , ignore_case = TRUE), '')
}


#
# Get features for a BigWig file
#
calc.features <- function(bw, type = 'max') {
  unlist(summary(bw, ranges.test.pp, type = type))
}
bw <- BigWigFile(bigwig.file)
gr <- ranges.test()[sample(length(ranges.test()), 10)]
feat <- unlist(summary(bw, gr, type = c('max')))
system.time(unlist(summary(bw, ranges.test.pp, type = c('max'))))


#
# Get the name of the feature file
#
features.file.name <- feature.file.name(feature)
if (file.exists(features.file.name)) {
  message('Features already exist, not recreating: ', features.file.name)
} else {
  #
  # Calculate and save features
  message('Calculating features without footprints')
  subsample <- function(gr) gr[sample(length(gr), 100000)]
  subsample(bw)
  length(bw)
  system.time(features <- lapply(lapply(lapply(bigwig, get.bw), subsample), calc.feature))
  # save.features(features, features.file.name)
}
message('Done')
