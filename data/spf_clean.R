##########################################################################
##### MSc Thesis - Niels Blom - ECB SPF Data Cleaning (Chapter 4)      ###
##########################################################################
#
# Stage 3 of the data pipeline (see README_data_pipeline.txt for Stages 1-2:
# merging raw SPF rounds via Lee & Seregina's Excel.py/Excel2.py, and
# updating realised GDP actuals from Eurostat).
#
# Input:  second_edit.csv                   (128 x 108, 1999Q3-2026Q2)
#         Actual_gdpgrowth_data_upd.xlsx    (107 quarters, 1999Q3-2026Q1)
# Output: yhat.csv, forERR.csv, ytrue.csv
#
##########################################################################

library(missForest)
library(readxl)
library(here)

###################################################################
##### PART 1: Load and Transpose                               #####
#####   second_edit.csv comes out of the Python merge step as  #####
#####   forecasters x rounds; transpose to the rounds x        #####
#####   forecasters (T x p) orientation used throughout the    #####
#####   empirical script                                       #####
###################################################################

raw    <- read.csv(here("data", "raw", "second_edit.csv"), row.names = 1, check.names = FALSE)
rounds <- colnames(raw)                  # extract round names before transposing
Y      <- t(raw)                         # forecasters x rounds -> rounds x forecasters
Y      <- apply(Y, 2, as.numeric)
rownames(Y) <- rounds                    # t() and as.numeric() both drop row names

cat("Raw: T =", nrow(Y), ", p =", ncol(Y), "\n")
cat("Coverage:", rounds[1], "to", rounds[length(rounds)], "\n")

###################################################################
##### PART 2: Missingness Filter                                #####
#####   Drop forecasters missing more than 40% of rounds -      #####
#####   irregular reporters would otherwise dominate the        #####
#####   imputation step below with very little real signal      #####
###################################################################

miss_frac <- colMeans(is.na(Y))
Y         <- Y[, miss_frac <= 0.40, drop = FALSE]
cat("After 40% filter: p =", ncol(Y), "\n")

###################################################################
##### PART 3: Random Forest Imputation                          #####
#####   missForest (Stekhoven & Buhlmann, 2012) fills the       #####
#####   remaining gaps for forecasters who report most, but not #####
#####   all, rounds - default settings used throughout          #####
###################################################################

set.seed(2025)
Y_imputed <- missForest(Y)$ximp
rownames(Y_imputed) <- rounds            # missForest drops row names; restore them
cat("Missing after imputation:", sum(is.na(Y_imputed)), "\n")

###################################################################
##### PART 4: Load Actuals and Align                            #####
#####   Match each survey round to its realised GDP growth      #####
#####   figure (Eurostat namq_10_gdp, EA20 - see README)        #####
###################################################################

actuals      <- read_excel(here("data", "raw", "Actual_gdpgrowth_data_upd.xlsx"))
ytrue_lookup <- setNames(actuals$y, actuals$`TIME PERIOD`)
ytrue        <- ytrue_lookup[rounds]

cat("Rounds without actuals:", sum(is.na(ytrue)),
    "->", rounds[is.na(ytrue)], "(will be dropped)\n")

###################################################################
##### PART 5: Drop Unmatched Rounds and Compute Forecast Errors #####
###################################################################

keep      <- !is.na(ytrue)
Y_imputed <- Y_imputed[keep, , drop = FALSE]
ytrue     <- ytrue[keep]
rounds    <- rounds[keep]

cat("Final: T =", nrow(Y_imputed), ", p =", ncol(Y_imputed), "\n")
cat("Final coverage:", rounds[1], "to", rounds[length(rounds)], "\n")

# Sign convention: error = actual - forecast (not forecast - actual).
forERR <- -sweep(Y_imputed / 100, 1, ytrue, FUN = "-")

cat("Mean forecast error:", round(mean(forERR), 4), "\n")
cat("Mean off-diagonal error correlation:",
    round(mean(cor(forERR)[lower.tri(cor(forERR))]), 3), "\n")

###################################################################
##### PART 6: Save                                              #####
#####   Outputs feed directly into empirical_spf_gdp.R          #####
###################################################################

write.csv(Y_imputed / 100, here("data", "yhat.csv"),  row.names = TRUE)
write.csv(forERR,          here("data", "forERR.csv"), row.names = TRUE)
write.csv(data.frame(round = rounds, actual = ytrue), here("data", "ytrue.csv"), row.names = FALSE)

cat("Done. Saved yhat.csv, forERR.csv, ytrue.csv\n")
