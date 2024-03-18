### R variables for sourcing - keep these
# dependencies
require("data.table")
require("stringr")

# base dirs
bd <- "~/Documents/MS/"
bda <- paste0(bd, "analysis/")
bds <- paste0(bd, "scripts/")
bdr <- paste0(bd, "raw/")
bdd <- paste0(bd, "doc/")
bdpdb <- paste0(bd, "pdb/")

# pipeline inputs --> anonymised
lcmsms89 <- paste0(bdr, "#########################5x/combined/txt/")
lcmsmsFAIMS89 <- paste0(bdr, "#########################5x_FAIMS/combined/txt/")
lcmsms99 <- paste0(bdr, "#########################4xNuclearLysate_MQ2-1-0-0/combined/txt/")

mq_pg_list_89 <- c(paste0(lcmsms89,"proteinGroups.txt"), paste0(lcmsmsFAIMS89,"proteinGroups.txt")) |> setattr(name = "filetype", value = "mq-pg-89")
#mq_pg_list_all <- c(mq_pg_list_89, paste0(lcmsms99,"proteinGroups.txt")) |> setattr(name = "filetype", value = "mq-pg-all")

mq_evi_list_89 <- c(paste0(lcmsms89,"evidence.txt"), paste0(lcmsmsFAIMS89,"evidence.txt")) |> setattr(name = "filetype", value = "mq-evi-89")
#mq_evi_list_all <- c(mq_pg_list_89, paste0(lcmsms99,"evidence.txt")) |> setattr(name = "filetype", value = "mq-evi-all")

mq_pg_list_99 <- paste0(lcmsms99,"proteinGroups.txt") |> setattr(name = "filetype", value = "mq-pg-99")
mq_evi_list_99 <- paste0(lcmsms99,"evidence.txt") |> setattr(name = "filetype", value = "mq-evi-99")


# output dirs
out_wide <- paste0(bda, "wide/")
    out_wide_evi <- paste0(out_wide, "wide/evidence/")
    out_wide_pg <- paste0(out_wide, "wide/proteinGroups/")
out_img <- paste0(bd, "img_volcano/")
out_volcano_pre <- paste0(bda, "volcano_pre/") 
out_reports <- paste0(bda, "reports/")

# for msgfplus
mzid_dir <-paste0(bda, "mzid/")
mzml_dir <-paste0(bda, "mzML/")
tsv_combi_dir <- paste0(bda, "mzid_datasets_combined/")


################################################################################
### Global settings - these can be modified
# data_preformatter_wide vars
#imputation_mode <- "QRILC" # QRILC (imputes for cols where at least 50% of the row is present), razor (discards all rows w/ NAs), null (does nothing besides ensuring that at least 1 value is present per sample)
peptide_length_min <- 6 # for manual protein-level intensity generation
peptide_length_max <- 25 # for manual protein-level intensity generation
prenormalisation_logbase <- 2
protein_choice <- "2" # 1 (protein with most peptides from protein groups) or 2 (a test; protein found most often in groups in the sample - this is unfinished and should throw an error if this is different from the default ordering); shouldn't really matter
#unique_peptides_min <- 3 # only proteins with at least this number of minimum unique peptides are kept
verbose <- F
writefile <- T

# volcano_preformatter vars
maxP <- 0.05
minFC <- 1.5 
cpus <- 8
#p_test <- "ROTS" # ROTS = Reproducibility Optimised Test Statistic, t-test = T-test with Benjamini-Hochberg correction



### for testing... #############################################################
# unique_peptides_min <- 3
# imputation_mode <- "razor"
# p_test <- "ROTS"
# session_name <- paste("PLK99",
#                       paste0(unique_peptides_min,"uniquePeptides"),
#                       paste0("log", prenormalisation_logbase),
#                       paste0("impute-", imputation_mode),
#                       paste0("test-", p_test),
#                       paste0("maxP", maxP),
#                       paste0("minFC", minFC),
#                       sep = "_")
# verbose <- F
# intensity_logbase <- 2
# p_test_FDR_cutoff <- 0.01
# ##############################################################################

# Render the Rmarkdown file to HTML; use different setting combinations and render the results
if (!exists("flag")) { # done to avoid cyclic rendering error; flag is set in the main code before sourcing this file
    testing <- "not testing" # pop testing flag; if not mentioned, beepr will notify at the end of each volcano pre/plot batch run
    force_overwrite <- F
    iteration_counter <- 0
    for (imputation_mode in c("QRILC", "razor", "duplicate")) {
        imputation_mode <- imputation_mode
        for (unique_peptides_min in c(2, 3)) {
            unique_peptides_min <- unique_peptides_min
            for (p_test in c("ROTS", "t-test")) {
                # # toy data
                # imputation_mode <- "QRILC"
                # unique_peptides_min <- 2
                # p_test <- "ROTS"
                
                # setup output vars
                session_name <- paste("PLK89",
                                      paste0(unique_peptides_min,"uniquePeptides"),
                                      paste0("log", prenormalisation_logbase),
                                      paste0("impute-", imputation_mode),
                                      paste0("test-", p_test),
                                      paste0("maxP", maxP),
                                      paste0("minFC", minFC),
                                      sep = "_")
                outdir_report <- paste0(out_reports, session_name, "/")
                out_file <- paste0(outdir_report, session_name, ".html")
                
                # prep output directory
                if (!dir.exists(outdir_report)) dir.create(outdir_report, recursive = T)
                
                # Check if the output file exists and its modification time
                if (file.exists(out_file)) {
                    file_mod_time <- file.info(out_file)$mtime
                    hours_diff <- as.numeric(difftime(Sys.time(), file_mod_time, units = "hours"))
                    
                    if (hours_diff <= 3 & !force_overwrite) {
                        # Skip rendering if the file was created in the last 3 hours
                        cat("\nSkipping rendering for", out_file, "as it was created less than 3 hours ago.\n")
                        iteration_counter <- iteration_counter + 1
                        next
                    }
                }
                
                # render file, update iterator
                cat("\n\nRendering file:", session_name)
                rmarkdown::render(input = paste0(bds, "MS_2023_binding_analysis.Rmd"), output_file = out_file)
                iteration_counter <- iteration_counter + 1
            }
        }
    }
    if (iteration_counter == 12) beepr::beep(4) else beepr::beep(6)
}