library(data.table)
setwd("~/julia_workflow")

reps <- 30

for (warming_level in c("0.0", "2.5", "5.0", "7.5", "10.0")) {
  for (i in 1:reps) {
    print(paste("Warming Level:", warming_level, "| Rep:", i))
    disp_model.df <- fread(paste("./outputs/disp_outputfile_", warming_level, "_C_", i, ".csv", sep = ""), header = TRUE, sep = ',')
    save(disp_model.df, file = paste("./outputs/disp_outputfile_", warming_level, "_C_", i, ".RData", sep =""))
  }
}

for (warming_level in c("0.0", "2.5", "5.0", "7.5", "10.0")) {
  for (i in 1:reps) {
    print(paste("Warming Level:", warming_level, "| Rep:", i))
    alpha_model.df <- fread(paste("./outputs/alpha_outputfile_", warming_level, "_C_", i, ".csv", sep = ""), header = TRUE, sep = ',')
    save(alpha_model.df, file = paste("./outputs/alpha_outputfile_", warming_level, "_C_", i, ".RData", sep =""))
  }
}

for (warming_level in c("0.0", "2.5", "5.0", "7.5", "10.0")) {
  for (i in 1:reps) {
    print(paste("Warming Level:", warming_level, "| Rep:", i))
    alpha_disp_model.df <- fread(paste("./outputs/alpha_disp_outputfile_", warming_level, "_C_", i, ".csv", sep = ""), header = TRUE, sep = ',')
    save(alpha_disp_model.df, file = paste("./outputs/alpha_disp_outputfile_", warming_level, "_C_", i, ".RData", sep =""))
  }
}

for (warming_level in c("0.0", "2.5", "5.0", "7.5", "10.0")) {
  for (i in 1:reps) {
    print(paste("Warming Level:", warming_level, "| Rep:", i))
    alpha_disp_model.df <- fread(paste("./outputs/asymmetry_outputfile_", warming_level, "_C_", i, ".csv", sep = ""), header = TRUE, sep = ',')
    save(alpha_disp_model.df, file = paste("./outputs/asymmetry_outputfile_", warming_level, "_C_", i, ".RData", sep =""))
  }
}



#files <- list.files("./outputs/", full.names = TRUE)

# Define a function to rename files recursively
#rename_files <- function(files) {
 # for (file in files) {
  #  if (tools::file_ext(file) == "csv") {
      # Rename the file
   #   new_file <- gsub("disp_", "alpha_", file)
    #  file.rename(file, new_file)
      #cat("Renamed", file, "to", new_file, "\n")
   # } else if (is.dir(file)) {
      # Recursively call the function for subdirectories
     # rename_files(list.files(file, full.names = TRUE))
   # }
 # }
#}

# Rename files recursively
#rename_files(files)
