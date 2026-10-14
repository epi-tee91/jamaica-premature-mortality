# Objective 1: Temporal Trends and Patterns

This folder contains the analysis plan and Stata code used to describe temporal trends and sociodemographic patterns in overall and cause-specific premature mortality in Jamaica between 2010 and 2019.

Premature mortality was defined as death before age 70 and measured primarily using years of potential life lost (YPLL).

## Contents

- `01_Indicators_Age70_Public.do`: Stata program used to prepare the analytic dataset and calculate mortality and YPLL indicators, directly age-standardized rates, standard errors, supplementary regression models, JCF homicide comparisons, and files used for Joinpoint regression.
- `Objective1_Analysis_Plan_Public.pdf`: Final dataset creation and analysis plan for Objective 1.

## Software

The primary analyses were conducted using Stata version 19. Temporal trend analyses were conducted using the National Cancer Institute’s Joinpoint Regression Program version 6.0.

## Running the program

Before running the Stata program, users must update the following line to identify the location of the project files on their computer:

```stata
global thesis_root "C:/path/to/Thesis"
