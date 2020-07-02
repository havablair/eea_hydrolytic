#### METADATA ####

# Purpose of script: process nitrate data for 2019 Exp1 PMN samples

# By: Miriam Gieske
# Last updated: 10/29/2019

#### LOAD PACKAGES ####
library(plyr)
library(tidyverse)
library(plater)
library(rio)

#### SET WORKING DIRECTORY ####
# working directory = folder where the files with data can be found
setwd("G:/HORT/GrossmanLab/UMN - Grossman Lab/Projects/2016 Grant Awards/OREI-HighTunnel-Anne, Julie/Data and R code/Exp1/PMN/Nitrate data and code 2019 and redos/Templates with data")
setwd("C:/Users/gies0107/Documents/OREI copy no photos 3-12-2020/Data and R code/Exp1/PMN/Nitrate data and code 2019 and redos/Templates with data")

#### CONVERT EXCEL FILES TO CSV ####

# You should only need to do this once (unless you add more excel files)

# Files needed:
# One completed template file (2019-08-16 version) for each plate, 
#   in xlsx format.
# Files used with this version of the script should contain all the metadata in the csv file itself,
#   in 5 blocks:
#     "Template" = sample and standard IDs.  Empty wells (no sample, standard, or check)
#          may be coded as ".", "0" (zero), "NA", or left blank.
#     "Data" = raw absorbance data from the plate reader
#     "Bad_wells" = identifies any wells that have known problems (ex. pipetting errors)
#     "Dilution" has the dilution factor for each sample
#     "Batch" identifies which samples are associated with which controls (must leave NAs blank)
# The files you want to process must all be in a single folder, 
#   and that folder must not contain any other xlsx files.  
#   R will attempt to read in all files in the designated folder.

# Convert the data to csv format so read.plates can read it in.
xls <- dir(pattern = "xlsx")

created <- mapply(convert, xls, gsub("xlsx", "csv", xls),
       MoreArgs = list(in_opts = list(which = "Template and data")))

#### READ IN THE CSV FILES ####

# If you get an error: check whether you have any extra csv files 
#   in the designated folder
file.names <- list.files(pattern = ".csv$", full.names = FALSE)
plates <- read_plates(file.names)
names(plates) <- c("Plate", "Wells", "ID", "abs", "Bad_wells", "Dilution", "Batch")
head(plates)
unique(plates$Plate) # To check that all plates read in correctly

#### REMOVE EMPTY WELLS AND KNOWN BAD DATA ####

# Remove empty wells  
plates <- subset(plates, !is.na(ID) & ID != 0 & ID !=".")

# Remove bad wells and keep good ones
plates <- subset(plates, Bad_wells != "bad" & Bad_wells != "Bad" 
                 & Bad_wells != "x" & Bad_wells != "X" | is.na(Bad_wells))

#### MEAN ABSORBANCE FOR EACH SAMPLE AND STANDARD ####

# Get mean of the three tech reps for each sample
mean.abs <- plates %>%
  group_by(Plate, ID) %>%
  summarise(abs = mean(abs, na.rm = T),
            dil = mean(Dilution),
            batch = mean(Batch))

head(mean.abs)

#### QUALITY CONTROL ####

# Count number of technical replicates for each sample (after bad wells are removed)
good.ct <- plates %>%
  group_by(Plate, ID) %>%
  summarise(good = sum(!is.na(abs)))

# Coefficient of variation (CV) is standard deviation divided by mean,
#   expressed as a percent
# Set the maximum acceptable CV for technical replicates of the same sample.
# Samples with a CV above this threshold will be removed from your dataset.
# Ex. 10 = remove samples with a CV greater than 10%
cv_max <- 10

# Calculate CV for each sample
cv.dat <- plates %>%
  group_by(Plate, ID) %>%
  summarise(sample.cv = sd(abs, na.rm = TRUE)/mean(abs, na.rm = TRUE)*100)

# Look at the distribution of CVs for each plate
# You can use this information to adjust your CV threshold, if desired
dlply(cv.dat, .(Plate), function(dat) 
  hist(dat$sample.cv, main = unique(dat$Plate), xlab = "CV"))

# Merge the quality control info with the mean absorbance dataset
mean.abs <- merge(mean.abs, good.ct, by = c("Plate", "ID"))
mean.abs <- merge(mean.abs, cv.dat, by = c("Plate", "ID"))
head(mean.abs)

# Add 2 columns to track QC data
mean.abs$keep <- "yes"   # Good data: yes, keep it.  Bad data: no, remove it.
mean.abs$reason <- NA    # Fill in for bad data.

# Check CV and flag samples that are over the threshold you set above
mean.abs$keep <- ifelse(mean.abs$sample.cv > cv_max, "no", mean.abs$keep)
mean.abs$reason <- ifelse(mean.abs$sample.cv > cv_max, "High CV", mean.abs$reason)

# Check number of good replicates and flag samples with less than 2
mean.abs$keep <- ifelse(mean.abs$good < 2, "no", mean.abs$keep)
mean.abs$reason <- ifelse(mean.abs$good < 2, "Pipetting errors", mean.abs$reason)

# Show the bad samples
mean.abs[mean.abs$keep=="no", ]

#### MAKE STANDARD CURVES ####

# Set the shared pattern of letters that identifies your standards
std_key <- "Std"

# Get the absorbance of the standards
stds.abs <- mean.abs[grepl(std_key, mean.abs$ID, fixed=TRUE)==TRUE 
                     & mean.abs$keep=="yes", ]

# Get the concentrations from the standard names
stds.abs$conc <- as.numeric(str_split_fixed(stds.abs$ID, "-", n=2)[ , 2])

# Remove standards that don't have a concentration value
stds.abs <- subset(stds.abs, !is.na(stds.abs$conc))

# Split into a list with one dataframe for each plate
stds.list <- split(stds.abs, stds.abs$Plate)

# Function to get equation and R-sq for a model
lm_eqn <- function(m) {
  l <- list(intcpt = format(abs(coef(m)[1]), digits = 3, scientific=FALSE),
            slope = format(coef(m)[2], digits = 3, scientific=FALSE),
            rsq = format(summary(m)["r.squared"], digits = 3));
  if (coef(m)[1] >= 0)  {
    eq <- paste("y = ", l["slope"], "x + ", l["intcpt"], ", R^2 = ", l["rsq"], sep = "")
  } else {
    eq <- paste("y = ", l["slope"], "x - ", l["intcpt"], ", R^2 = ", l["rsq"], sep = "")    
  }
}

# Function to plot a standard curve with equation and R-sq
plot.curve <- function(data, model) {
  plot(abs ~ conc, data, main = unique(data$plate)) +
  abline(model) +
  mtext(lm_eqn(model), side = 3)
}

# Set up matrix to record coefficients and max conc for each plate's standard curve
coeff <- matrix(data = NA, ncol = 4, nrow = length(stds.list))
colnames(coeff)=c("Plate", "intcpt", "slope", "max.abs")

# Set up a matrix to record notes for each plate
curve.notes <- matrix(data = NA, ncol = 2, nrow = length(stds.list))
colnames(curve.notes)=c("Plate", "notes")

# Set the minimum number of standards for a curve
# The absolute minimum to plot a straight line is 3
# We normally use 5-6 standards
min_stds <- 4

# This loop does the following:
# 1. Omits plates with fewer than the minimum number of standards passing QC checks
# 2. Plots standard curves with equation and R-sq for each remaining plate
# 3. Returns coefficients for plates with R-sq > 0.99
# 4. Prints outcome for each plate to the console

for(i in 1:length(stds.list)) {
  # Check number of standards
  if(length(stds.list[[i]]$ID) < min_stds) {
    # Record the plate name
    coeff[i, "Plate"] <- names(stds.list)[i]
    # Feedback
    curve.notes[i, ] <- c(names(stds.list[i]), "Need to redo plate (not enough good standards)")
  } else {
    # Calculate curve
    m1 <- lm(abs ~ conc, stds.list[[i]])
    rsq1 <- summary(m1)["r.squared"] 
    # Plot the standard curve
    plot.curve(stds.list[[i]], m1)
    if (rsq1 >= 0.98) {
      # Save the coefficients
      coeff[i, "Plate"] <- names(stds.list)[i]
      coeff[i, 2:3] <- coef(m1)
      coeff[i, 4] <- max(stds.list[[i]]$abs)
      # Feedback
      curve.notes[i, ] <- c(names(stds.list[i]), "Standard curve looks good!")
    } else {
        # Record the plate name
        coeff[i, "Plate"] <- names(stds.list)[i]
        # Feedback
        curve.notes[i, ] <- c(names(stds.list[i]), "Need to redo plate (bad standard curve)")
      }
      # Clean up
    rm(m1)
    rm(rsq1)
  }
}

print(curve.notes)

#### ADD INFO TO THE SAMPLE DATASET ####

# Remove standards from absorbance dataset
samples <- mean.abs[grepl(std_key, mean.abs$ID, fixed=TRUE)==FALSE, ]

# Add coefficients and max conc for standard curves
coeff <- data.frame(coeff, stringsAsFactors = FALSE)
samples <- merge(samples, coeff, by = "Plate")
samples$intcpt <- as.numeric(samples$intcpt)
samples$slope <- as.numeric(samples$slope)
samples$max.abs <- as.numeric(samples$max.abs)

# Flag data from plates with bad standard curves
samples$keep <- ifelse(is.na(samples$intcpt), "no", samples$keep)
samples$reason <- ifelse(is.na(samples$intcpt), "Bad standard curve", samples$reason)

# If the absorbance of a sample is greater than the absorbance
#   of the highest concentration standard, flag the data
samples$keep <- ifelse(!is.na(samples$max.abs) & samples$abs > samples$max.abs, 
                       "no", samples$keep)
samples$reason <- ifelse(!is.na(samples$max.abs) & samples$abs > samples$max.abs, 
                         "Outside linear range", samples$reason)

#### LOOK AT THE ABSORBANCE OF YOUR CHECKS ####

# How does the absorbance of the soil-less checks compare to the absorbance
#   of the 0 ppm standards?
# Note: I do not want to include the plates I know to be bad in this analysis

# Extract the 0 ppm standards
zero.stds <- stds.abs[stds.abs$conc == 0 & stds.abs$keep == "yes",
                      c("ID", "abs", "keep", "reason")]
zero.stds$type <- "0ppm_std"

# Set the shared pattern of letters that allows you to identify checks
check_key <- "CK"
# Extract checks from dataset
checks <- samples[grepl(check_key, samples$ID, fixed=TRUE)==TRUE 
                  & samples$keep == "yes", 
                  c("ID", "abs", "keep", "reason")]
checks$type <- "soil-less_check"

# Make histograms
# If your checks are clean, the absorbance distribution for the checks should look
#   similar to the absorbance distribution for the 0ppm standards.
stds.checks <- rbind(zero.stds, checks)
ggplot(stds.checks, aes(abs, fill = type)) + 
  geom_histogram(alpha = 0.5, position = 'identity', binwidth = 0.001) +
  theme_classic() +
  ggtitle("Nitrate 2019")

# If you have checks with high absorbance, take a look at them 
#   to see which ones are high.  High absorbance values for checks 
#   suggest contamination.  You should carefully evaluate whether the data
#   you collected along with those checks is usable.
checks2 <- samples[grepl(check_key, samples$ID, fixed=TRUE)==TRUE, 
                   c("Plate", "ID", "batch", "abs", "keep", "reason")]
checks2[order(checks2$abs),] 

#### CALCULATE CONCENTRATION OF SAMPLES ####
# For POX-C, this will be the conc. of permanganate remaining unreduced after the reaction.
# Note: in the equation for the standard curve, x = conc and y = abs
#   abs = slope*conc + intcpt
# We have abs and need to solve for conc
#   conc = (abs - intcpt)/slope
# We also need to take into account the dilution factor
samples$conc <- with(samples, ((abs - intcpt)/slope)*dil)

#### SAVE THE RESULTS ####

# Set the folder where you want to save the results
setwd("G:/HORT/GrossmanLab/UMN - Grossman Lab/Projects/2016 Grant Awards/OREI-HighTunnel-Anne, Julie/Data and R code/Exp1/PMN/Nitrate data and code 2019 and redos")

# Save the good data
good.samples <- subset(samples, keep=="yes")
write.csv(good.samples, "OREI_Exp1_nitrate_2019_processed.csv", row.names = F)

# Save the samples that need to be redone
bad.samples <- subset(samples, keep=="no")
write.csv(bad.samples, "OREI_Exp1_nitrate_2019_bad.csv", row.names = F)
