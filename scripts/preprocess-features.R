#!/usr/bin/env Rscript
#
# Load RDS features objects and preprocess with as(..., 'Matrix') then save again
#
#!
#! sbatch directives begin here ###############################
#!
#SBATCH -p mrc-bsu-sand                 # Partition
#SBATCH -J PREPFEAT                     # Name of the job
#SBATCH -A MRC-BSU-SL2                  # Which project should be charged
#SBATCH --nodes=1                       # How many whole nodes should be allocated?
#SBATCH --ntasks=1                      # How many (MPI) tasks will there be in total? (<= nodes*16)
#SBATCH --mem=16000                     # How many MB each node is allocated
#SBATCH --time=06:00:00                 # How much wallclock time will be required?
#SBATCH -o prepfeat-%j.out              # stdout
#SBATCH -e prepfeat-%j.out              # stderr
#SBATCH --mail-type=FAIL                # What types of email messages do you wish to receive?
##SBATCH --no-requeue                   # Uncomment this to prevent the job from being requeued (e.g. if
"Usage: preprocess-features.R FEAT..." -> doc

#
# Set warnings as errors
#
options(warn = 2)


#
# Load packages
#
devtools::load_all()


#
# Parse options
#
# .args <- "../Data/Motifs/Known/scan-gr.rds"
if (! exists(".args")) .args <- commandArgs(TRUE)
opts <- docopt::docopt(doc, args = .args)
print(opts)
rds.files <- opts[['FEAT']]


#
# Convert the objects if necessary
#
for (rds.file in rds.files) {
  message('Loading   : ', rds.file)
  obj <- readRDS(rds.file)
  if (! is(obj, 'Matrix')) {
    message('Converting: ', rds.file)
    obj.pp <- as(obj, 'Matrix')
    message('Saving    : ', rds.file)
    saveRDS(obj.pp, rds.file)
  }
}
message('Done')
