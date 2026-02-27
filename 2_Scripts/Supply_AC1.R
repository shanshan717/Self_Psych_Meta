# Environment setup
base::rm(list = base::ls()) 
base::library(tidyverse)
base::library(irrCAC)

# Define file paths
input_path  <- "1_Data/AnalysisData/Literature_info/Interrater_Reliability_Raw.csv"
output_path <- "3_Output/7_other/Interrater_Reliability_Processed.csv"

# Data Import and Preprocessing

# Import the raw reliability data
# check.names = FALSE ensures column names like 'A1', 'A2' are preserved
df_irr <- readr::read_csv(input_path, col_names = TRUE, show_col_types = FALSE) %>% 
  base::as.data.frame()

# Filter out empty rows often found at the end of CSV files
df_irr <- df_irr %>%
  dplyr::filter(!base::is.na(Variable) & Variable != "")

# Automated Consistency Scoring (0/1)

# Identify columns starting with 'A' followed by numbers (e.g., A1, A2...)
all_colnames <- base::names(df_irr)
a_cols <- base::grep("^A[0-9]+$", all_colnames, value = TRUE)

# Sort columns numerically to ensure correct pairing (A1 with A2, A3 with A4)
a_cols <- a_cols[base::order(base::as.numeric(base::sub("A", "", a_cols)))]

# Loop through pairs to calculate binary agreement scores
r_count <- 1
for (i in base::seq(1, base::length(a_cols), by = 2)) {
  
  # Define pair columns
  col_coder1 <- a_cols[i]
  col_coder2 <- a_cols[i + 1]
  r_col_name <- base::paste0("R", r_count)
  
  # Logical comparison: 1 if identical (case-insensitive), 0 otherwise
  df_irr[[r_col_name]] <- base::ifelse(
    base::trimws(base::tolower(base::as.character(df_irr[[col_coder1]]))) == 
    base::trimws(base::tolower(base::as.character(df_irr[[col_coder2]]))),
    1, 0
  )
  
  # Handle potential NA values by assigning them as 0 (disagreement)
  df_irr[[r_col_name]][base::is.na(df_irr[[r_col_name]])] <- 0
  
  r_count <- r_count + 1
}


# ====================================================================
#  Inter-rater Reliability Analysis (AC1)
# ====================================================================

# Select only the generated agreement columns (R1, R2, R3...)
ratings_matrix <- df_irr %>% 
  dplyr::select(dplyr::matches("^R\\d+"))

base::print(base::paste("Number of agreement items extracted:", base::ncol(ratings_matrix)))

# Flatten matrix to a vector for Gwet's AC1 calculation
# Following Sun et al. (2022), compare actual scores against a perfect standard (all 1s)
actual_vector <- base::as.vector(base::as.matrix(ratings_matrix))
actual_vector <- actual_vector[!base::is.na(actual_vector)]

cac_input <- base::data.frame(
  Standard = base::rep(1, base::length(actual_vector)), 
  Actual   = actual_vector
)

# Calculate AC1 using irrCAC package
ac1_result <- irrCAC::gwet.ac1.raw(cac_input)

# Results Reporting

base::cat("\n========================================\n")
base::cat(" INTER-RATER RELIABILITY (AC1) REPORT \n")
base::cat("========================================\n")

# Safely extract coefficients and confidence intervals
coeff <- base::as.numeric(ac1_result$est$coeff.val)
lcb   <- base::as.numeric(ac1_result$est$conf.int[1])
ucb   <- base::as.numeric(ac1_result$est$conf.int[2])
pa    <- base::as.numeric(ac1_result$est$pa)

base::cat("AC1 Coefficient:      ", base::round(coeff, 4), "\n")
base::cat("95% Conf. Interval:   [", base::round(lcb, 4), ",", base::round(ucb, 4), "]\n")
base::cat("Overall Agreement (%): ", base::round(pa, 4) * 100, "%\n")
base::cat("========================================\n")

# Detailed result table check
# base::print(ac1_result$est)