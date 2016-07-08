################################################################################
## Download competition data for the
## ENCODE-DREAM in vivo Transcription Factor Binding Site Prediction Challenge
################################################################################
library(synapseClient)   # Synapse Client, 1.13-4

## If you haven't set up a hook in your .synapseConfig file, you'll have to
## login to Synapse:
#synapseLogin(username = 'me@example.com', password = 'secret')

cat("Make sure you've accepted the terms of use before running this script!\n")

## You may wish to copy these files to a specific destination directory. If so,
## set the path to that directory here:
dest_dir = 'Data'

# ChIPseq fold_change_signal = syn6181334
# ChIPseq labels = syn6181335
# ChIPseq peaks = syn6181336
# DNASE bams = syn6176232
# DNASE fold_coverage_wiggles = syn6176233
# DNASE peaks = syn6176234
# RNAseq = syn6176231
folder_ids = c('syn6181334', 'syn6181335', 'syn6181336', 'syn6176232', 
               'syn6176233', 'syn6176234', 'syn6176231')

data_files = list()

for (folder_id in folder_ids) {
     
     ## get folder
     folder <- synGet(folder_id)
     cat(sprintf('Downloading contents of %s folder (%s)\n',
                 propertyValue(folder, 'name'), propertyValue(folder, 'id')))
     
     ## query for child entities
     query_results = synapseQuery(sprintf('select id,name from file where parentId=="%s"', 
                                          folder_id))
     
     for (i in seq(nrow(query_results))) {
          name <- query_results[i,'file.name']
          id <- query_results[i, 'file.id']
          
          ## get file extension
          parts <- strsplit(name, '\\.')[[1]]
          extension <- parts[length(parts)]
          
          cat(sprintf('\tDownloading file: %s (%s)\n', name, id))
          data_files <- append(data_files, synGet(id))
          
     }
}

cat('\ndownload complete!\n')

if (!is.null(dest_dir)) {
     cat(sprintf('Copying files to %s\n', dest_dir))
     for (file in data_files) {
          cat(sprintf('\tCopying %s\n', basename(getFileLocation(file))))
          file.copy(getFileLocation(file), dest_dir)
     }
}

cat('\ncopying complete!\n')
