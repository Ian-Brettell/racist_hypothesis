#!/usr/bin/env Rscript


#############################
# Libraries
#############################

library(here)
library(tidyverse)
library(pegas)
library(knitr)
library(ggridges)
library(gwasrapidd)
library(cowplot)
library(plotly)
library(DT)

#############################
# Paths
#############################

## Latest plot path
plot_path = here::here("plots", "20210203_batch")
dir.create(plot_path, showWarnings = F)

#############################
# Functions
#############################

## Read `.afreq` files from `Plink`
read_afreq <- function(file){
  out = read.table(file, header = T, comment.char = "") %>%
    dplyr::rename(CHR = X.CHROM,
                  SNP = ID) %>% 
    dplyr::mutate(CHR = as.integer(CHR))

  return(out)
}

#############################
# Parameters
#############################

is.color = function (x) 
{
  if (is.null(x)) 
    return(FALSE)
  return(setNames(vapply(x, function(X) {
    tryCatch(is.matrix(grDevices::col2rgb(X)), error = function(e) FALSE)
  }, FUN.VALUE = TRUE), NULL))
}

darker = function (col, amount = 150) 
{
  if (!all(is.color(col))) 
    stop("All elements in col must be valid colors. Use is.col(col) to check it.")
  if (!methods::is(amount, "numeric") || length(amount) != 
      1) 
    stop("amount must be a single number")
  if (amount > 255 || amount < 0) 
    stop("amount must be a number between 0 and 255")
  .darker <- function(col, amount) {
    if (is.na(col)) 
      return(NA)
    new.col <- ((grDevices::col2rgb(col)) - amount)/255
    new.col[new.col[, 1] < 0, 1] <- 0
    return(grDevices::rgb(t(new.col)))
  }
  return(unlist(lapply(col, .darker, amount)))
}

lighter = function (col, amount = 150) 
{
  if (!all(is.color(col))) 
    stop("All elements in col must be valid colors. Use is.col(col) to check it.")
  if (!methods::is(amount, "numeric") || length(amount) != 
      1) 
    stop("amount must be a single number")
  if (amount > 255 || amount < 0) 
    stop("amount must be a number between 0 and 255")
  .lighter <- function(col, amount) {
    if (is.na(col)) 
      return(NA)
    new.col <- ((grDevices::col2rgb(col)) + amount)/255
    new.col[new.col[, 1] > 1, 1] <- 1
    return(grDevices::rgb(t(new.col)))
  }
  return(unlist(lapply(col, .lighter, amount)))
}

target_traits = c("body height",
                  "body mass index",
                  "self reported educational attainment",
                  "intelligence",
                  "inflammatory bowel disease",
                  "skin pigmentation",
                  "skin pigmentation measurement",
                  "eye color",
                  "eye colour measurement",
                  "hair color",
                  "hair colour measurement",
                  "schizophrenia",
                  "unipolar depression",
                  "fasting blood glucose measurement",
                  "myocardial infarction",
                  "low density lipoprotein cholesterol measurement",
                  "platelet count")


# Factor levels for `HIT_CONTROL`
hit_control_levels = c("hit", "control")

# Create vectors for recoding traits with full names and vice versa
recode_vec = c("hei" = "Height",
               "bmi" = "BMI",
               "edu" = "Educational attainment",
               "int" = "Intelligence",
               "ibd" = "IBD",
               "pig" = "Pigmentation")

rev_recode_vec = names(recode_vec)
names(rev_recode_vec) = recode_vec

# Create vectors to recode pigmentation traits
pig_traits = c("skin pigmentation",
               "skin pigmentation measurement",
               "eye color",
               "eye colour measurement",
               "hair color",
               "hair colour measurement")
pig_recode_vec = rep("all pigmentation", length(pig_traits))
names(pig_recode_vec) = pig_traits

# Colour palettes
pal_primary = c("Height" = "#FC4E07",
                "BMI" = "#FFBF00",
                "Educational attainment" = "#0BC166",
                "Intelligence" = "#00AFBB",
                "IBD" = "#D84797",
                "Pigmentation" = "#360568")
pal_secondary = c("Height" = "#B63502",
                  "BMI" = "#B88A00",
                  "Educational attainment" = "#07743D",
                  "Intelligence" = "#00727A",
                  "IBD" = "#972062",
                  "Pigmentation" = "#1F033A")

onekg_pal = c("#B8984D","#E3B64C","#CFB54D","#D49943","#C26B36","#E1B759","#ECCB51","#682B7B","#8C5793","#8D3B5E", "#AE307F","#5D448A","#B8C650","#7FAA53","#8DAF4F","#5E8A48","#6E974B","#2D3468","#394D92","#798EC1","#95C4DB","#81B6C2","#B6302C","#B1253A","#A33E3A","#A34028")
names(onekg_pal) = c("LWK", "GWD", "MSL", "ACB", "ASW", "YRI", "ESN", "BEB", "STU", "ITU", "PJL", "GIH", "CHB", "KHV", "CHS", "JPT", "CDX", "TSI", "CEU", "IBS", "GBR", "FIN", "PEL", "MXL", "CLM", "PUR")


# New pal with extended traits
extended_traits = c("body height",
                    "body mass index",
                    "self reported educational attainment",
                    "intelligence",
                    "inflammatory bowel disease",
                    "all pigmentation",
                    "schizophrenia",
                    "unipolar depression",
                    "fasting blood glucose measurement",
                    "myocardial infarction",
                    "low density lipoprotein cholesterol measurement",
                    "platelet count" )

pal_primary_new = c("#fc4e07","#ffbf00","#0bc166","#00afbb","#d84797","#360568",
                    "#4f0943", "#c200fb", "#00647a", "#57C13A", "#f2a918", "#e84141")

names(pal_primary_new) = extended_traits

# New pal with different shade for each pigmentation trait

pal_primary_new_exp = c("#fc4e07","#ffbf00","#0bc166","#00afbb","#d84797",
                      lighter("#360568", amount = 80),
                      lighter("#360568", amount = 60),
                      lighter("#360568", amount = 40),
                      lighter("#360568", amount = 20),
                      "#360568",
                      darker("#360568", amount = 20),
                      "#4f0943", "#c200fb", "#00647a", "#57C13A", "#f2a918", "#e84141")
names(pal_primary_new_exp) = target_traits

# Pal for 
traits_with_pig = c("body height",
                    "body mass index",
                    "self reported educational attainment",
                    "intelligence",
                    "inflammatory bowel disease",
                    "skin pigmentation / measurement",
                    "schizophrenia",
                    "unipolar depression",
                    "fasting blood glucose measurement",
                    "myocardial infarction",
                    "low density lipoprotein cholesterol measurement",
                    "platelet count" )
pal_primary_with_pig = pal_primary_new
names(pal_primary_with_pig) = traits_with_pig

# Turbo palette

source("https://gist.githubusercontent.com/jlmelville/be981e2f36485d8ef9616aef60fd52ab/raw/466a6564a86066a9600860cf8058fab5d23e1da5/turbo_colormap.R")
