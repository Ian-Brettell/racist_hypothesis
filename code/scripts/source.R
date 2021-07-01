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

## Get manhattan plot
get_man <- function(df, trait, title, chr, bp, snp, p){
  # Plot
  manhattan(df, chr=chr, bp=bp, snp=snp, p=p,
            cex = 0.6,
            c(pal_primary[names(pal_primary) == trait],
              pal_secondary[names(pal_secondary) == trait]),
            main = title)
}


#############################
# Parameters
#############################

# Import karyoploteR `lighter` and `darker` functions

source("https://gist.githubusercontent.com/brettellebi/c5015ee666cdf8d9f7e25fa3c8063c99/raw/91e601f82da6c614b4983d8afc4ef399fa58ed4b/karyoploteR_lighter_darker.R")


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

## Clumping target parameter
clump_param = "r2-0.1_kb-1000"

# Factor levels for `trait` (or `PHENO`)
trait_levels = c("hei", "bmi", "edu", "int", "ibd", "pig")
names(trait_levels) = trait_levels

trait_levels_verb =  c("Height",
                       "BMI",
                       "Educational attainment",
                       "Intelligence",
                       "IBD",
                       "Pigmentation")
names(trait_levels_verb) = trait_levels_verb

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
