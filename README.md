# PinealGlandChicken

This repository was created aiming to make all the data available to ensure the reproducibility of the scientific paper entitled:

Methylation, miRNA and gene expression have an integrated sex-specific role in the pineal gland of birds subjected to unpredictable light schedules

And published in:
(www. ???)

RawData Availability: 
The dataset supporting the conclusions of this article is available from the European Nucleotide Archive (ENA) repository (EMBL-EBI), under accession number PRKEB35831 (www.ebi.ac.uk/ena/data/view/PRKEB35831).

From the count tables generated from the raw data for each omic level analyzed, the scripts available here allow the reproducibility of this study.

They were divided into three parts:

Part 1: Statistical Analysis for the 5 contrasts used in each omic level (DMR, DMiR, DEG) using Limma (https://bioconductor.org/packages/release/bioc/html/limma.html)
Contrasts:
comp1= M_Stress-F_Stress ,
comp2= M_Control-F_Control ,
comp3= M_Stress-M_Control ,
comp4= F_Stress-F_Control ,

Foldchanges:
comp1, for exemple:
+ means = F>M
- means = F<M


2nd part: Scripts for loading all the inputs generated from 1nd part and counting features from different dataframes for tables generation.

3rd part: Scripts to make graphs and tables from the inputs of the 2nd part.
