# Expected workflow using CranioMorph
.rs.restartR(clean = TRUE)

pak::pak("/home/kristen/Documents/CranioMorph_0.0.1.9000.tar.gz", ask = FALSE)

library(CranioMorph)
rm(list = ls())

# Load data from servers - this will be different when our servers are reorganized
data = query_data(species = "mm", dataset = "Bmp")
data$meta = data$meta[data$meta$Gene == "Bmp7",]

# Create a "borealis" object from data
# These are default parameters. Adjust and provide to argument if you want to change it.
analysis_params = list(
  use_PC_scores = FALSE, # Use PCs to reduce shape dimensions?
  remove_size = FALSE, # remove the effects of csize on shape in data processing?
  assess_allometry = TRUE, # assess the effects of csize on shape in analyses?
  log_size = FALSE # log-transform csize?
)

borealis = initialise_borealis(data, analysis_params = analysis_params)

# TODO: interactive process to set up analysis_params (in initialise?)

# Specify columns if you don't want to go through the derive_meta interactions.
# IF EFFECTS OF SEX AND AGE SHOULD NOT BE REMOVED AS CONFOUNDERS, INPUT AS independant VARIABLES.

column_definitions = list(
  ID = 1,
  Sex = 5,
  Age = numeric(0),
  Confounder = numeric(0),
  independant = c(4)
)

# Detect metadata; use automated = FALSE to have specific control
borealis = derive_meta(borealis, automated = TRUE, column_definitions = column_definitions, impute = FALSE)
  # meta will now have identifiers as ROWNAMES
  # borealis object also contains $identifiers now

#TODO: if automated = TRUE and column_definitions provide - do not enter detecting logic and just push through the provided column definitions

# # After metadata is detected, subsetting using borealis objects by [individuals, landmarks]
# test = data[data$identifiers[1:50],1:3]

borealis = process_shapes(borealis, symmetrize = "detect", correction = "Experimental_Group", force = TRUE)

res = find_diffs(borealis, include_allometry = TRUE)
res$shapes

# look at final model
summary(borealis$analysis_params$correction_model$Control)

# saveRDS(borealis$shapes, "~/Downloads/testborealis.rds")
#
# #######
# library(dplyr)
#
# test = borealis[borealis$meta$independant_Gene %in% c("Bmp7", "Twsg1")]
# covs = borealis$meta %>%
#   select(sex, age)



# group =
#
#
# test = residualise(gpa$rotated, covs = covs, interactions = TRUE, group = borealis$meta$independant_Gene)
# print(test)
#
#
# shapes = gpa$rotated
# interactions = TRUE
# group = borealis$meta$independant_Gene

# sympair_748 = as.matrix(readRDS("/mnt/Storage4/Users/Kristen/bhlab_misc/landmark_sym.Rds")$sym.pairings)
# usethis::use_data(sympair_748, overwrite = TRUE)
