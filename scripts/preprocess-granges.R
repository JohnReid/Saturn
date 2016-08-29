#!/usr/bin/env Rscript
#
# Load RDS GRanges objects and preprocess with GNCList then save again
#
#!
#! sbatch directives begin here ###############################
#!
#SBATCH -p mrc-bsu-sand                 # Partition
#SBATCH -J PREPROCESS                   # Name of the job
#SBATCH -A MRC-BSU-SL2                  # Which project should be charged
#SBATCH --nodes=1                       # How many whole nodes should be allocated?
#SBATCH --ntasks=1                      # How many (MPI) tasks will there be in total? (<= nodes*16)
#SBATCH --mem=15000                     # How many MB each node is allocated
#SBATCH --time=192:00:00                # How much wallclock time will be required?
#SBATCH -o preprocess-%j.out            # stdout
#SBATCH -e preprocess-%j.out            # stderr
#SBATCH --mail-type=FAIL                # What types of email messages do you wish to receive?
##SBATCH --no-requeue                   # Uncomment this to prevent the job from being requeued (e.g. if
"Usage: preprocess-granges.R GR..." -> doc

#
# Set warnings as errors
#
options(warn = 2)


#
# Load packages
#
library(GenomicRanges)


#
# Parse options
#
# .args <- "../Data/Motifs/Known/scan-gr.rds"
if (! exists(".args")) .args <- commandArgs(TRUE)
opts <- docopt::docopt(doc, args = .args)
print(opts)
rds.files <- opts[['GR']]


#' Preprocess the argument if it is not already a GNCList
#'
asGNCList <- function(gr) {
  if (is(gr, "GNCList")) {
    gr
  } else {
    GNCList(gr)
  }
}


#
# Convert the ranges
#
for (rds.file in rds.files) {
  message('Loading: ', rds.file)
  gr <- readRDS(rds.file)
  if (is(gr, "list")) {
    gr.pp <- lapply(gr, asGNCList)
  } else if (is(gr, "GNCList")) {
    gr.pp <- asGNCList(gr)
  }
  message('Saving: ', rds.file)
  saveRDS(gr.pp, rds.file)
}
