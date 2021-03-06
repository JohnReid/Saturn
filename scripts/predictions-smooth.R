#!/usr/bin/env Rscript
#!
#! sbatch directives begin here ###############################
#!
#SBATCH -p mrc-bsu-sand                    # Partition
#SBATCH -J PREDSMOOTH                      # Name of the job
#SBATCH -A MRC-BSU-SL2                     # Which project should be charged
#SBATCH --nodes=1                          # How many whole nodes should be allocated?
#SBATCH --ntasks=1                         # How many (MPI) tasks will there be in total? (<= nodes*16)
#SBATCH --cpus-per-task=1                  # How many CPUs will be used per task
#SBATCH --mem=30400                        # How many MB each node is allocated
#SBATCH --time=06:00:00                    # How much wallclock time will be required?
#SBATCH -o pred-smooth/pred-smooth-%j.out  # stdout
#SBATCH -e pred-smooth/pred-smooth-%j.out  # stderr
#SBATCH --mail-type=FAIL                   # What types of email messages do you wish to receive?
##SBATCH --no-requeue                      # Uncomment this to prevent the job from being requeued (e.g. if
                                           # interrupted by node failure or system downtime):
"Usage: predictions-smooth.R [options] IN OUT

Options:
  --override=TF                  Use parameters for TF defined in R package [default: ""]
  -l --length-scale=LENGTHSCALE  Use LENGTHSCALE [default: 50]
  --logodds                      Smooth log-odds instead of probabilities [default: FALSE]
  --width=MAXWIDTH               Width of kernel in units of regions [default: 20]" -> doc


#
# Load package
#
devtools::load_all()
library(Saturn)


#
# Set warnings as errors
#
options(warn = 2)


#
# Parse options
#
# .args <- "--logodds ../Data/Predictions/predictions.xgboost.chrfold.EGR1.H1-hESC.DNase_Known_KnownWell_DREME_DREMEWell.tsv ../slurm/smoothed-predictions.tsv"
# .args <- '--length-scale=10 --width=5 ../Data/Predictions/predictions.xgboost.remove.FOXA2.HepG2.DNase_Known_KnownWell_DREME_DREMEWell.tsv ../Data/Predictions/predictions.xgboost.L=10-LO=no-W=5.FOXA2.HepG2.DNase_Known_KnownWell_DREME_DREMEWell.tsv'
# options(error = recover)
if (! exists(".args")) .args <- commandArgs(TRUE)  # Check if we have manually set arguments for debugging
opts <- docopt::docopt(doc, args = .args)
print(opts)
in.path <- opts$IN
out.path <- opts$OUT
override <- opts[['override']]
log.transform <- as.logical(opts[['logodds']])
max.width <- as.integer(opts[['width']])
length.scale <- as.numeric(opts[['length-scale']])


#
# Configure parameters
#
if ("" != override) {
  message('Overriding with parameters for TF: ', override)
  .params <- filter(smoothing.params, TF == override)
  stopifnot(1 == length(.params))
  log.transform <- as.logical(.params$LO)
  length.scale <- as.numeric(.params$L)
  max.width <- as.integer(.params$W)
  actually.smooth <- as.logical(.params$smooth)
} else {
  actually.smooth <- TRUE
}


#
# Are we actually going to do any smoothing?
#
if (! actually.smooth) {
  message('No need to smooth, will just copy predictions unchanged')
  file.copy(in.path, out.path)
  quit(save = "no")
}


#
# Show parameters
#
message('Log transform: ', log.transform)
message('Length scale: ', length.scale)
message('Maximum width: ', max.width)


#
# Load predictions
#
message('Loading predictions from: ', in.path)
preds <- data.table::fread(
  in.path,
  col.names = c('chrom', 'start', 'end', 'prediction', 'bound'))
sample_n(preds, 15)
sapply(preds, class)


#
# Smooth predictions
#
predictions.smoothed <- smooth.predictions(preds, length.scale, max.width, log.transform)


#
# Write predictions
#
message('Writing predictions to: ', out.path)
preds$prediction <- as.vector(predictions.smoothed)
data.table::fwrite(preds, out.path, sep = '\t', col.names = FALSE)
