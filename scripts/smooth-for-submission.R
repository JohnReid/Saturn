#!/usr/bin/env Rscript
#!
#! sbatch directives begin here ###############################
#!
#SBATCH -p mrc-bsu-sand                    # Partition
#SBATCH -J SUBSMOOTH                       # Name of the job
#SBATCH -A MRC-BSU-SL2                     # Which project should be charged
#SBATCH --nodes=1                          # How many whole nodes should be allocated?
#SBATCH --ntasks=1                         # How many (MPI) tasks will there be in total? (<= nodes*16)
#SBATCH --cpus-per-task=1                  # How many CPUs will be used per task
#SBATCH --mem=60800                        # How many MB each node is allocated
#SBATCH --time=06:00:00                    # How much wallclock time will be required?
#SBATCH -o sub-smooth/sub-smooth-%j.out    # stdout
#SBATCH -e sub-smooth/sub-smooth-%j.out    # stderr
#SBATCH --mail-type=FAIL                   # What types of email messages do you wish to receive?
##SBATCH --no-requeue                      # Uncomment this to prevent the job from being requeued (e.g. if
                                           # interrupted by node failure or system downtime):
"Usage: smooth-for-submission.R IN..." -> doc


#
# Load package
#
devtools::load_all()
library(Saturn)
library(stringr)


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


#
# For each input file
#
for (in.path in opts$IN) {
  message('Handling: ', in.path)
  #
  # Deduce TF and cell type from filename
  #
  tf.pred <- as.vector(str_match(in.path, tf.regex))
  cell.pred <- as.vector(str_match(in.path, cell.regex))
  message('TF   : ', tf.pred)
  message('Cell : ', cell.pred)
  tfs.tf.cell <- filter(tfs, TF == tf.pred, cell == cell.pred)
  print(tfs.tf.cell)


  #
  # Set output name
  #
  if ('submit' == tfs.tf.cell$split) {
    label <- 'F'
  } else if ('ladder' == tfs.tf.cell$split) {
    label <- 'L'
  } else {
    stop('Unknown split type.')
  }
  out.path = str_c(saturn.data(), '/Submissions/', label, '.', tf.pred, '.', cell.pred, '.tab')
  message('Output: ', out.path)


  #
  # Configure parameters
  #
  .params <- filter(smoothing.params, TF == tf.pred)
  stopifnot(1 == nrow(.params))
  log.transform <- as.logical(.params$LO)
  length.scale <- as.numeric(.params$L)
  max.width <- as.integer(.params$W)
  actually.smooth <- as.logical(.params$smooth)
  message('Log transform: ', log.transform)
  message('Length scale: ', length.scale)
  message('Maximum width: ', max.width)


  #
  # Load predictions
  #
  message('Loading predictions from: ', in.path)
  preds <- data.table::fread(in.path, header = FALSE)
  col.names = c('chrom', 'start', 'end', 'prediction', 'bound')
  colnames(preds) <- col.names[1:ncol(preds)]
  # sample_n(preds, 15)


  #
  # Are we actually going to do any smoothing?
  #
  if (! actually.smooth) {
    message('No need to smooth, will just copy predictions unchanged')
    predictions.smoothed <- preds$prediction
  } else {
    predictions.smoothed <- smooth.predictions(preds, length.scale, max.width, log.transform)
  }


  #
  # Write predictions
  #
  message('Writing predictions to: ', out.path)
  preds$prediction <- as.vector(predictions.smoothed)
  data.table::fwrite(preds, out.path, sep = '\t', col.names = FALSE)


  #
  # Compress file
  #
  message('Compressing predictions: ', out.path)
  system2('/bin/gzip', args = out.path)
}
