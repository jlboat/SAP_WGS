#!/usr/bin/env bash

R -e """
library(mashr)
library(data.table)

# Read and filter data
setwd('/zfs/tillers/lucas/SAP/gwas/QUALITY/recal')
df <- fread(file='PH.MLM.csv', header=T)
df <- df[!is.infinite(df$SE),]

# Generate ash data
ash_data <- ash(df$Effect, df$SE)

# Filter effects and SE by lfsr < 0.1
effects <- ash_data$data$x[get_lfsr(ash_data) < 0.1]
standard_error <- ash_data$data$s[get_lfsr(ash_data) < 0.1]

# Generate Mash Data
data = mash_set_data(effects, standard_error)

# Randomly Sample 100,000 markers for Mash control

"""
