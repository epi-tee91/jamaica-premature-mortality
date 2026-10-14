/*********************************************************************************************************************
File Name: 01_Indicators_Age70

Purpose: To produce the indicators for study objective 1 (**WITH A REFERENCE AGE OF 70 FOR DEFINING PREMATURE MORTALITY**). The indicators are:
			1. Cumulative Incidence of Mortality  
			2. Percentage of Deaths that are Premature     
			3. Years of Potential Life Lost (YPLL)
			4. Mean Years of Potential Life Lost (YPLL) per Death    
			5. Unconditional Probability of Dying Between Ages 30 And 70 Years from Non-Communicable Diseases
	
Author: Tonya Campbell

Begun: 2024-05-28

Completed: 2024-10-03

Revised: 2024-11-27 (added lines of code to extract datasets with standardized rates for joinpoint regression)
		 2025-04-29 (revising calculations based on updated CODs)
		 2025-05-22 (adding standard errors)
		 2025-08-04 (adding chi^2 tests, calculation of median ypll [see note starting on line 1857], and wilcoxon rank sum tests)
		 2026-05-14 (updating standard errors for standardized rates and adding for crude and age-specific rates)
		 2026-06-24 (adding supplementary table comparing JCF and RGD homicides)
		 2026-06-25 (revised SE calculation to Poisson approximation)
		2026-08-25 (added Poisson regression models to examine sex-based differences in population rates of YPLL and premature mortality, on recommendation from Prof YC)
		2026-09-17 (revised analysis comparing JCF and RGD homicides, added another analysis of JCF homicide trends)
		2026-09-21 (revised crude and directly age-standardized YPLL standard errors to use the sum of squared individual YPLL contributions.)
		2026-09-23 (added Table 1 of descriptive characteristics)

Statistical Reviews: 2024-11-25 by Novie Younger-Coleman
					2026-08-04 by Novie Younger-Coleman
					2026-09-21 — Final AI-assisted static epidemiologic and statistical code review: The complete Objective 1 age-70 indicator do-file was reviewed following revision of the 
								YPLL standard-error calculations. The review confirmed that crude and directly age-standardized YPLL variances now use the sum of squared individual
								YPLL contributions, consistent with a weighted-Poisson framework. Mortality, YPLL, standardization, NCD probability, JCF trend, hypothesis-testing, Joinpoint-input,
								and supplementary Poisson-regression sections were reviewed. Minor revisions were identified concerning explicit year and sex restrictions, variable-name 
								abbreviations, merge validation, terminology, sparse chi-square tests, and documentation of the expanded homicide definition.

Dataset creation plan: See the Objective 1 analysis plan included in the repository.

Input Datasets: deaths_2010_2020 (dataset of registered deaths from 2010 to 2020)
				population_overall (dataset holding values for the Jamaican population from 2010 to 2019)
				population_sex (dataset holding values for the Jamaican population by sex from 2010 to 2019)
				population_age (dataset holding values for the Jamaican population by age group from 2010 to 2019)
				population_agesex (dataset holding values for the Jamaican population by age group and sex from 2010 to 2019)
				population_age_5yrbands (dataset holding values for the Jamaican population by age group in 5-year bands from 2010 to 2019)
				population_agesex_5yrbands (dataset holding values for the Jamaican population by age group in 5-year bands and sex from 2010 to 2019)
				whostandard_all (dataset holding the WHO standard population weights for the entire population)
				whostandard_70 (dataset holding the WHO standard population weights for the population of people <70 years of age)
				whostandard_65 (dataset holding the WHO standard population weights for the population of people <65 years of age)
				jcf_murders_2013_2019.xlsx (dataset count of homicides recorded by JCF from 2013-2019 by age)

Output: 
		Datasets: obj1_deaths_2010_2019
				  obj1_overall70 
				  obj1_overall70_stdize 
				  obj1_sex70 
				  obj1_sex70_stdize
				  obj1_age70
				  obj1_agesex70
				  obj1_overallsex70_crude_jnpt
				  obj1_age70_crude_jnpt
				  obj1_overallsex70_stdize_jnpt
				  obj1_chi2results_70 
				  obj1_wmwresults_70
				  
				  
		Excel: Output - Age 70.xlsx
*********************************************************************************************************************/

version 19
clear all
set more off
capture log close _all

/*******************************************************************************
USER-DEFINED DIRECTORIES

Set thesis_root to the location of the main folder on the computer
running this program.

The individual-level mortality data are restricted and are not included in
the public repository.
*******************************************************************************/

global thesis_root "C:/path/to/Thesis" /* Users must replace the path with their local directory. */

global obj1_root             "${thesis_root}/Objective 1"
global clean_data_dir        "${thesis_root}/Data/03_Data Cleaning/04_Cleaned Data"
global analysis_dir          "${obj1_root}/03_Data Analysis"
global project_data_dir      "${analysis_dir}/Project Data"
global input_dir             "${analysis_dir}/Input"
global analysis_output_dir   "${analysis_dir}/Output"
global manuscript_output_dir "${obj1_root}/04_Output"

/*******************************************************************************
Set working directory, save log file
*******************************************************************************/
cd "${project_data_dir}"
log using "${analysis_dir}/01_Indicators_Age70.log", replace name(Log)



/*******************************************************************************
Load dataset
*******************************************************************************/
use "${clean_data_dir}/deaths_2010_2020.dta" 
describe 


/*******************************************************************************
Prepare Dataset
*******************************************************************************/


/*********Calculate YPLL*********/

/*reference age 70*/
generate ypll_70 = 70 - age if premature70==1
	replace ypll_70 = . if premature70==0

*check for errors
assert !missing(ypll_70) if premature70 == 1
assert missing(ypll_70) if premature70 == 0

summarize ypll_70

*create a variable that replaces missing ypll with 0 (for use in standardized rate calculations)
gen ypll_70_burden = ypll_70
replace ypll_70_burden = 0 if missing(ypll_70_burden)

*create squared individual YPLL contribution for variance calculations
gen double ypll_70_sq_burden = ypll_70_burden^2

label variable ypll_70_sq_burden "Squared individual YPLL contribution, reference age 70"

assert ypll_70_sq_burden == 0 if premature70 == 0
assert ypll_70_sq_burden == ypll_70^2 if premature70 == 1
assert ypll_70_sq_burden >= 0


/*********COMBINE UNDERLYING CAUSE OF DEATH INDICATOR VARIABLES *(EXCEPT FOR THE COMPOSITE 4-CAUSE NCD VARIABLE)* INTO A SINGLE VARIABLE*********/

/*quality check -  # of causes of death per person - no one should have >1 cod*/
egen codcheck=rowtotal(underlying_cause_matperi underlying_cause_infect underlying_cause_cancer underlying_cause_diabendo underlying_cause_cvd underlying_cause_resp underlying_cause_digest underlying_cause_gu underlying_cause_psyc underlying_cause_external underlying_cause_other)

tab1 codcheck, missing  /*looks okay - proceed with categorizing*/

gen underlying_cause_cat = .
local i = 1
foreach v of varlist underlying_cause_matperi underlying_cause_infect underlying_cause_cancer underlying_cause_diabendo underlying_cause_cvd underlying_cause_resp underlying_cause_digest underlying_cause_gu underlying_cause_psyc underlying_cause_external underlying_cause_other {
    replace underlying_cause_cat = `i' if `v' == 1
    label define underlying_cause_cat `i' "`v'", add
    local ++i /*i+1*/
}

label define underlying_cause_catlbl 1 "Maternal & Perinatal Conditions" 2 "Infectious & Parasitic Diseases" 3 "Cancer" 4 "Diabetes & Endocrine Disorders" 5 "Circulatory Diseases" 6 "Respiratory Diseases" 7 "Digestive Diseases" 8 "Genitourinary Diseases" 9 "Neuropsychiatric Conditions" 10 "External Causes" 11 "Other Causes", replace
label values underlying_cause_cat underlying_cause_catlbl

tab1 underlying_cause_cat, missing


**restrict to 2010-2019**
keep if inrange(death_year, 2010, 2019)
assert inrange(death_year, 2010, 2019)

tab1 death_year, missing


/*********Other minor prep*********/

/*create a temporary variable with all values = 1 to allow for counting total number of deaths*/
generate _temp = 1 

rename death_year year


/*********Apply Exclusions*********/

**exclude people who resided outside of Jamaica at the time of death***
tab1 _temp, missing
tab1 residence_parish, missing
count if !inrange(residence_parish,1,14)
keep if inrange(residence_parish,1,14)
tab1 residence_parish, missing nolabel 

count if missing(age)
drop if missing(age)


save obj1_deaths_2010_2019, replace

sort year sex underlying_cause_cat
by year: summarize ypll_70, detail
by year sex underlying_cause_cat: summarize ypll_70, detail
clear

/*******************************************************************************
Indicators to reported according to the following groups (check DCP for 
additional details):
	- Overall
	- By Age Group
	- By Sex
	- By Age Group for each Sex
	
For each indicator, results will be reported for all causes and by cause of death
*******************************************************************************/


/*******************************************************************************
REFERENCE AGE 70
*******************************************************************************/

/************************************
Get stats for entire study period
************************************/
use obj1_deaths_2010_2019

/* Standardize ICD-10 codes */
capture drop cod_check r_code

gen str10 cod_check = upper( ///
    subinstr( ///
    subinstr(trim(underlying_cause_clean), ".", "", .), ///
    " ", "", .))

/* Identify R00–R99 */
gen byte r_code = regexm(cod_check, "^R[0-9][0-9]")

/* Count and percentage in final analytic sample */
quietly count if r_code == 1
local n_r = r(N)

quietly count
local n_total = r(N)

display "R-code deaths: " `n_r'
display "Total deaths: " `n_total'
display "Percentage R-code: " %5.2f (100 * `n_r' / `n_total')

/* Confirm assigned cause category */
tab underlying_cause_cat if r_code == 1, missing

use obj1_deaths_2010_2019, clear

collapse (count) total_mortality_count=_temp (mean) premature_percentage=premature70 (sum) premature_count=premature70 total_ypll_count=ypll_70 (median) median_ypll=ypll_70 (mean) mean_ypll=ypll_70



/************************ 
 SECTION 1 - OVERALL
************************/

use obj1_deaths_2010_2019, clear

tab1 age, missing

/****SECTION 1.1 - OVERALL POPULATION, ALL CAUSE MORTALITY - CRUDE RATES****/

/*use collapse function to get summary statistics, which will be the numerators for indicators #1-4...
		indicator #5 will be programmed later*/

collapse /// 
	(count) total_mortality_count = _temp ///
	(mean) premature_percentage = premature70 ///
	(sum) premature_count = premature70 /// 
	(sum) total_ypll_count = ypll_70 ///
	(sum) total_ypll_squared = ypll_70_sq_burden ///
	(median) median_ypll = ypll_70 /// 
	(mean) mean_ypll = ypll_70, ///
 by(year)   

replace premature_percentage = premature_percentage * 100

format total_mortality_count total_ypll_count premature_count %9.0gc
format premature_percentage %8.2f
format median_ypll mean_ypll %4.2f

/*Add in population data and calculate rates*/

*overall population

frame rename default overall

frame create populationdata
frame change populationdata

use "${clean_data_dir}/population_overall.dta" 

frame change overall

frlink 1:1 year, frame(populationdata)

frget population, from(populationdata)

*population under age 70

frame create populationdata70
frame change populationdata70

use "${clean_data_dir}/population_agesex.dta" 

tab1 agegroup

drop if inlist(agegroup,7,8)

tab1 agegroup

collapse (sum) population70=population, by(year) 

frame change overall

frlink 1:1 year, frame(populationdata70)

frget population70, from(populationdata70)

*calculate rates of mortality and premature mortality 

generate total_mortality_rate = (total_mortality_count/population)*100000
generate total_ypll_rate = (total_ypll_count/population70)*100000
generate se_ypll_rate = (sqrt(total_ypll_squared)/population70)*100000

drop total_ypll_squared

format total_mortality_rate total_ypll_rate se_ypll_rate %9.0gc

order year population population70 total_mortality_count total_mortality_rate premature_count premature_percentage total_ypll_count total_ypll_rate se_ypll_rate median_ypll mean_ypll

rename (total_mortality_count total_mortality_rate premature_count premature_percentage total_ypll_count total_ypll_rate se_ypll_rate median_ypll mean_ypll) (total_mortality_count_all total_mortality_rate_all premature_count_all premature_percentage_all total_ypll_count_all total_ypll_rate_all se_ypll_rate_all median_ypll_all mean_ypll_all)

drop populationdata populationdata70

/*indicator 5 - unconditional probability of death between ages 30 and 70 from NCDs*/
frame create ncdprob
frame change ncdprob

use obj1_deaths_2010_2019

collapse (sum) total_mortality_count_ncd=underlying_cause_sdgncd, by(year agegroup_5yr)

tab1 agegroup, missing
tab1 agegroup, nolabel 

drop if agegroup==. | agegroup<7 | agegroup>14 /*probability to be calculated between ages 30 and 70*/

tab1 agegroup, missing

/*merge in population data*/
frame create population5yrbands
frame change population5yrbands

use "${clean_data_dir}/population_age_5yrbands.dta" 

frame change ncdprob

frlink 1:1 year agegroup_5yr, frame(population5yrbands)

frget population, from(population5yrbands)

drop population5yrbands

/*calculate unconditional probability of dying*/
generate mortrate = total_mortality_count_ncd / population /*step 1: calculate age-specific mortality rates for each five-year age group*/

generate mortprob = mortrate*5 / (1+mortrate*2.5) /*step 2: translate the mortality rate into the probability of death in each 5-yer age range*/

/*step 3: multiply the probabilities across age groups to get the total probability*/
keep agegroup_5yr year mortprob total_mortality_count_ncd

reshape wide total_mortality_count_ncd mortprob, i(year) j(agegroup_5yr)

generate probability_death_ncd = 1-((1-mortprob7)*(1-mortprob8)*(1-mortprob9)*(1-mortprob10)*(1-mortprob11)*(1-mortprob12)*(1-mortprob13)*(1-mortprob14))

/*obtain the count to correspond with the percentage*/
generate mortality_count_ncd = total_mortality_count_ncd7 + total_mortality_count_ncd8 + total_mortality_count_ncd9 + total_mortality_count_ncd10 + total_mortality_count_ncd11 + total_mortality_count_ncd12 + total_mortality_count_ncd13 + total_mortality_count_ncd14

format probability_death_ncd %04.3f
format mortality_count_ncd %9.0gc

/*link back to original frame*/

frame change overall

frlink 1:1 year, frame(ncdprob)

frget mortality_count_ncd probability_death_ncd, from(ncdprob)

drop ncdprob


/****SECTION 1.2 - OVERALL POPULATION, ALL CAUSE MORTALITY - AGE-STANDARDIZED RATES****/

/*Total Mortality*/

/*get counts*/
frame create stdizemortall
frame change stdizemortall

use obj1_deaths_2010_2019

collapse (count) total_mortality_count=_temp, by(year agegroup_5yr)   

drop if agegroup==.

/*merge in population data*/
frlink 1:1 year agegroup_5yr, frame(population5yrbands)

frget population, from(population5yrbands)

drop population5yrbands

/*merge in WHO standard population weights*/
frame create weightsall
frame change weightsall

use "${clean_data_dir}/whostandard_all.dta" 

frame change stdizemortall

frlink m:1 agegroup_5yr, frame(weightsall)

frget weight_whostd_all, from(weightsall)

drop weightsall

/*calculate standardized rates*/
generate crude = total_mortality_count / population

generate product = crude*weight_whostd_all
gen var_component = (weight_whostd_all^2) * (total_mortality_count / (population^2))
 
egen _total_mortality_rate_adj_all = total(product), by(year)
egen _se_adj_all = total(var_component), by(year)

generate total_mortality_rate_adj_all = _total_mortality_rate_adj_all*100000
generate se_adj_all=sqrt(_se_adj_all) * 100000

contract year total_mortality_rate_adj_all se_adj_all

format total_mortality_rate_adj_all %9.0gc


/*Premature Mortality*/

/*get counts*/
frame create stdizeypllall
frame change stdizeypllall

use obj1_deaths_2010_2019

collapse ///
    (sum) total_ypll_count = ypll_70_burden ///
	(sum) total_ypll_squared = ypll_70_sq_burden, ///
    by(year agegroup_5yr)

drop if agegroup_5yr==.

/*merge in population data*/

frlink 1:1 year agegroup_5yr, frame(population5yrbands)

frget population, from(population5yrbands)

drop population5yrbands

/*merge in WHO standard population weights (for subset of pop <70 years)*/
frame create weights70
frame change weights70

use "${clean_data_dir}/whostandard_70.dta" 

frame change stdizeypllall

frlink m:1 agegroup_5yr, frame(weights70)

frget weight_whostd_70, from(weights70)

drop weights70

tab1 agegroup_5yr, nolabel 
drop if agegroup_5yr >= 15 

/*calculate standardized rates*/
generate crude = total_ypll_count / population
generate product = crude*weight_whostd_70

gen double var_component = (weight_whostd_70^2 * total_ypll_squared) / (population^2)

/*sum across age groups*/
collapse ///
    (sum) _total_ypll_rate_adj_all = product ///
    (sum) _se_ypll_adj_all = var_component, ///
    by(year)

generate total_ypll_rate_adj_all = _total_ypll_rate_adj_all * 100000
gen se_ypll_adj_all = (sqrt(_se_ypll_adj_all)) * 100000

format total_ypll_rate_adj_all se_ypll_adj_all %9.0gc



/****SECTION 1.3 - OVERALL POPULATION STRATIFIED BY CAUSE OF DEATH - CRUDE RATES****/

frame create overallbymjrcause
frame change overallbymjrcause

use obj1_deaths_2010_2019

collapse /// 
	(count) total_mortality_count = _temp ///
	 (mean) premature_percentage = premature70 ///
	 (sum) premature_count = premature70 ///
	 (sum) total_ypll_count = ypll_70 ///
	 (sum) total_ypll_squared = ypll_70_sq_burden ///
	 (median) median_ypll = ypll_70 ///
	 (mean) mean_ypll = ypll_70, ///
	 by(year underlying_cause_cat) 

replace premature_percentage = premature_percentage * 100

format total_mortality_count total_ypll_count premature_count %9.0gc
format premature_percentage %8.2f
format median_ypll mean_ypll %4.2f

drop if underlying_cause_cat==.

reshape wide total_mortality_count total_ypll_count total_ypll_squared premature_count premature_percentage median_ypll mean_ypll, i(year) j(underlying_cause_cat)

/*Add in population data and calculate rates*/

*overall population

frlink 1:1 year, frame(populationdata)
frget population, from(populationdata)

*population under age 70

frlink 1:1 year, frame(populationdata70)
frget population70, from(populationdata70)

*calculate rates of mortality and premature mortality, rename and label variables

local codcat "matperi infect cancer diabendo cvd resp digest gu psyc external other"
local j = 1
forvalues j = 1/11 {
	local sffx : word `j' of `codcat' 
	rename (total_mortality_count`j' premature_percentage`j' premature_count`j' total_ypll_count`j' total_ypll_squared`j' median_ypll`j' mean_ypll`j') (total_mortality_count_`sffx' premature_percentage_`sffx' premature_count_`sffx' total_ypll_count_`sffx' total_ypll_squared_`sffx' median_ypll_`sffx' mean_ypll_`sffx')
    local ++j
}

local codcat2 "matperi infect cancer diabendo cvd resp digest gu psyc external other"
local q = 1
foreach x of local codcat2 {
    generate total_mortality_rate_`x' = (total_mortality_count_`x'/population)*100000
	generate total_ypll_rate_`x' = (total_ypll_count_`x'/population70)*100000
	generate se_ypll_rate_`x' = (sqrt(total_ypll_squared_`x')/population70)*100000
	format total_mortality_rate_`x' total_ypll_rate_`x' se_ypll_rate_`x' %9.0gc
    local ++q
}

drop total_ypll_squared_*

/*link frames with the crude rates, label variables*/

frlink 1:1 year, frame(overall)

frget total_mortality_count_all total_mortality_rate_all premature_count_all premature_percentage_all total_ypll_count_all total_ypll_rate_all se_ypll_rate_all median_ypll_all mean_ypll_all mortality_count_ncd probability_death_ncd, from(overall)

label variable population70 "Jamaica Mid-Year Population Estimates (Pop <70 Years of Age)"

label variable total_mortality_count_all "Total Mortality Count - Overall"
label variable total_mortality_rate_all "Total Mortality Rate - Overall"
label variable premature_count_all "No. Premature Deaths - Overall"
label variable premature_percentage_all "Pctg Premature Deaths - Overall"
label variable total_ypll_count_all "Total YPLL Count - Overall"
label variable total_ypll_rate_all "Total YPLL Rate - Overall"
label variable median_ypll_all "Median YPLL/Person - Overall"
label variable mean_ypll_all "Mean YPLL/Person - Overall"
label variable mortality_count_ncd "Number of Premature NCD Deaths (WHO Defn)"
label variable probability_death_ncd "Probability of Premature NCD Mortality (WHO Defn)"

label variable total_mortality_count_matperi "Total Mortality Count - Maternal & Perinatal Conditions"
label variable total_mortality_rate_matperi "Total Mortality Rate - Maternal & Perinatal Conditions"
label variable premature_count_matperi "No. Premature Deaths - Maternal & Perinatal Conditions"
label variable premature_percentage_matperi "Pctg Premature Deaths - Maternal & Perinatal Conditions"
label variable total_ypll_count_matperi "Total YPLL Count - Maternal & Perinatal Conditions"
label variable total_ypll_rate_matperi "Total YPLL Rate - Maternal & Perinatal Conditions"
label variable median_ypll_matperi "Median YPLL/Person - Maternal & Perinatal Conditions"
label variable mean_ypll_matperi "Mean YPLL/Person - Maternal & Perinatal Conditions"

label variable total_mortality_count_infect "Total Mortality Count - Infectious Diseases"
label variable total_mortality_rate_infect "Total Mortality Rate - Infectious Diseases"
label variable premature_count_infect "No. Premature Deaths - Infectious Diseases"
label variable premature_percentage_infect "Pctg Premature Deaths - Infectious Diseases"
label variable total_ypll_count_infect "Total YPLL Count - Infectious Diseases"
label variable total_ypll_rate_infect "Total YPLL Rate - Infectious Diseases"
label variable median_ypll_infect "Median YPLL/Person - Infectious Diseases"
label variable mean_ypll_infect "Mean YPLL/Person - Infectious Diseases"

label variable total_mortality_count_cancer "Total Mortality Count - Cancer"
label variable total_mortality_rate_cancer "Total Mortality Rate - Cancer"
label variable premature_count_cancer "No. Premature Deaths - Cancer"
label variable premature_percentage_cancer "Pctg Premature Deaths - Cancer"
label variable total_ypll_count_cancer "Total YPLL Count - Cancer"
label variable total_ypll_rate_cancer "Total YPLL Rate - Cancer"
label variable median_ypll_cancer "Median YPLL/Person - Cancer"
label variable mean_ypll_cancer "Mean YPLL/Person - Cancer"

label variable total_mortality_count_diabendo "Total Mortality Count - Diabetes & Diabetes & Endocrine Disorders"
label variable total_mortality_rate_diabendo "Total Mortality Rate - Diabetes & Endocrine Disorders"
label variable premature_count_diabendo "No. Premature Deaths - Diabetes & Endocrine Disorders"
label variable premature_percentage_diabendo "Pctg Premature Deaths - Diabetes & Endocrine Disorders"
label variable total_ypll_count_diabendo "Total YPLL Count - Diabetes & Endocrine Disorders"
label variable total_ypll_rate_diabendo "Total YPLL Rate - Diabetes & Endocrine Disorders"
label variable median_ypll_diabendo "Median YPLL/Person - Diabetes & Endocrine Disorders"
label variable mean_ypll_diabendo "Mean YPLL/Person - Diabetes & Endocrine Disorders"

label variable total_mortality_count_cvd "Total Mortality Count - Circulatory Diseases"
label variable total_mortality_rate_cvd "Total Mortality Rate - Circulatory Diseases"
label variable premature_count_cvd "No. Premature Deaths - Circulatory Diseases"
label variable premature_percentage_cvd "Pctg Premature Deaths - Circulatory Diseases"
label variable total_ypll_count_cvd "Total YPLL Count - Circulatory Diseases"
label variable total_ypll_rate_cvd "Total YPLL Rate - Circulatory Diseases"
label variable median_ypll_cvd "Median YPLL/Person - Circulatory Diseases"
label variable mean_ypll_cvd "Mean YPLL/Person - Circulatory Diseases"

label variable total_mortality_count_resp "Total Mortality Count - Respiratory Diseases"
label variable total_mortality_rate_resp "Total Mortality Rate - Respiratory Diseases"
label variable premature_count_resp "No. Premature Deaths - Respiratory Diseases"
label variable premature_percentage_resp "Pctg Premature Deaths - Respiratory Diseases"
label variable total_ypll_count_resp "Total YPLL Count - Respiratory Diseases"
label variable total_ypll_rate_resp "Total YPLL Rate - Respiratory Diseases"
label variable median_ypll_resp "Median YPLL/Person - Respiratory Diseases"
label variable mean_ypll_resp "Mean YPLL/Person - Respiratory Diseases"

label variable total_mortality_count_digest "Total Mortality Count - Digestive Diseases"
label variable total_mortality_rate_digest "Total Mortality Rate - Digestive Diseases"
label variable premature_count_digest "No. Premature Deaths - Digestive Diseases"
label variable premature_percentage_digest "Pctg Premature Deaths - Digestive Diseases"
label variable total_ypll_count_digest "Total YPLL Count - Digestive Diseases"
label variable total_ypll_rate_digest "Total YPLL Rate - Digestive Diseases"
label variable median_ypll_digest "Median YPLL/Person - Digestive Diseases"
label variable mean_ypll_digest "Mean YPLL/Person - Digestive Diseases"

label variable total_mortality_count_gu "Total Mortality Count - Genitourinary Diseases"
label variable total_mortality_rate_gu "Total Mortality Rate - Genitourinary Diseases"
label variable premature_count_gu "No. Premature Deaths - Genitourinary Diseases"
label variable premature_percentage_gu "Pctg Premature Deaths - Genitourinary Diseases"
label variable total_ypll_count_gu "Total YPLL Count - Genitourinary Diseases"
label variable total_ypll_rate_gu "Total YPLL Rate - Genitourinary Diseases"
label variable median_ypll_gu "Median YPLL/Person - Genitourinary Diseases"
label variable mean_ypll_gu "Mean YPLL/Person - Genitourinary Diseases"

label variable total_mortality_count_psyc "Total Mortality Count - Neuropsychiatric Conditions"
label variable total_mortality_rate_psyc "Total Mortality Rate - Neuropsychiatric Conditions"
label variable premature_count_psyc "No. Premature Deaths - Neuropsychiatric Conditions"
label variable premature_percentage_psyc "Pctg Premature Deaths - Neuropsychiatric Conditions"
label variable total_ypll_count_psyc "Total YPLL Count - Neuropsychiatric Conditions"
label variable total_ypll_rate_psyc "Total YPLL Rate - Neuropsychiatric Conditions"
label variable median_ypll_psyc "Median YPLL/Person - Neuropsychiatric Conditions"
label variable mean_ypll_psyc "Mean YPLL/Person - Neuropsychiatric Conditions"

label variable total_mortality_count_external "Total Mortality Count - External Causes"
label variable total_mortality_rate_external "Total Mortality Rate - External Causes"
label variable premature_count_external "No. Premature Deaths - External Causes"
label variable premature_percentage_external "Pctg Premature Deaths - External Causes"
label variable total_ypll_count_external "Total YPLL Count - External Causes"
label variable total_ypll_rate_external "Total YPLL Rate - External Causes"
label variable median_ypll_external "Median YPLL/Person - External Causes"
label variable mean_ypll_external "Mean YPLL/Person - External Causes"

label variable total_mortality_count_other "Total Mortality Count - Other"
label variable total_mortality_rate_other "Total Mortality Rate - Other"
label variable premature_count_other "No. Premature Deaths - Other"
label variable premature_percentage_other "Pctg Premature Deaths - Other"
label variable total_ypll_count_other "Total YPLL Count - Other"
label variable total_ypll_rate_other "Total YPLL Rate - Other"
label variable median_ypll_other "Median YPLL/Person - Other"
label variable mean_ypll_other "Mean YPLL/Person - Other"

label variable se_ypll_rate_all "Standard Error for Crude YPLL Rate - Overall"
label variable se_ypll_rate_matperi "Standard Error for Crude YPLL Rate - Maternal & Perinatal Conditions"
label variable se_ypll_rate_infect "Standard Error for Crude YPLL Rate - Infectious Diseases"
label variable se_ypll_rate_cancer "Standard Error for Crude YPLL Rate - Cancer"
label variable se_ypll_rate_diabendo "Standard Error for Crude YPLL Rate - Diabetes & Endocrine Disorders"
label variable se_ypll_rate_cvd "Standard Error for Crude YPLL Rate - Circulatory Diseases"
label variable se_ypll_rate_resp "Standard Error for Crude YPLL Rate - Respiratory Diseases"
label variable se_ypll_rate_digest "Standard Error for Crude YPLL Rate - Digestive Diseases"
label variable se_ypll_rate_gu "Standard Error for Crude YPLL Rate - Genitourinary Diseases"
label variable se_ypll_rate_psyc "Standard Error for Crude YPLL Rate - Neuropsychiatric Conditions"
label variable se_ypll_rate_external "Standard Error for Crude YPLL Rate - External Causes"
label variable se_ypll_rate_other "Standard Error for Crude YPLL Rate - Other"


order year population population70 total_mortality_count_all total_mortality_rate_all total_mortality_count_matperi total_mortality_rate_matperi total_mortality_count_infect total_mortality_rate_infect total_mortality_count_cancer total_mortality_rate_cancer total_mortality_count_diabendo total_mortality_rate_diabendo total_mortality_count_cvd total_mortality_rate_cvd total_mortality_count_resp total_mortality_rate_resp total_mortality_count_digest total_mortality_rate_digest total_mortality_count_gu total_mortality_rate_gu total_mortality_count_psyc total_mortality_rate_psyc total_mortality_count_external total_mortality_rate_external total_mortality_count_other total_mortality_rate_other premature_count_all premature_percentage_all premature_count_matperi premature_percentage_matperi premature_count_infect premature_percentage_infect premature_count_cancer premature_percentage_cancer premature_count_diabendo premature_percentage_diabendo premature_count_cvd premature_percentage_cvd premature_count_resp premature_percentage_resp premature_count_digest premature_percentage_digest premature_count_gu premature_percentage_gu premature_count_psyc premature_percentage_psyc premature_count_external premature_percentage_external premature_count_other premature_percentage_other total_ypll_count_all total_ypll_rate_all total_ypll_count_matperi total_ypll_rate_matperi total_ypll_count_infect total_ypll_rate_infect total_ypll_count_cancer total_ypll_rate_cancer total_ypll_count_diabendo total_ypll_rate_diabendo total_ypll_count_cvd total_ypll_rate_cvd total_ypll_count_resp total_ypll_rate_resp total_ypll_count_digest total_ypll_rate_digest total_ypll_count_gu total_ypll_rate_gu total_ypll_count_psyc total_ypll_rate_psyc total_ypll_count_external total_ypll_rate_external total_ypll_count_other total_ypll_rate_other median_ypll_all median_ypll_matperi median_ypll_infect median_ypll_cancer median_ypll_diabendo median_ypll_cvd median_ypll_resp median_ypll_digest median_ypll_gu median_ypll_psyc median_ypll_external median_ypll_other mean_ypll_all mean_ypll_matperi mean_ypll_infect mean_ypll_cancer mean_ypll_diabendo mean_ypll_cvd mean_ypll_resp mean_ypll_digest mean_ypll_gu mean_ypll_psyc mean_ypll_external mean_ypll_other mortality_count_ncd probability_death_ncd se_ypll_rate_all se_ypll_rate_matperi se_ypll_rate_infect se_ypll_rate_cancer se_ypll_rate_diabendo se_ypll_rate_cvd se_ypll_rate_resp se_ypll_rate_digest se_ypll_rate_gu se_ypll_rate_psyc se_ypll_rate_external se_ypll_rate_other

drop populationdata populationdata70 overall

save obj1_overall70, replace


/****SECTION 1.4 - OVERALL POPULATION, STRATIFIED BY CAUSE OF DEATH - AGE-STANDARDIZED RATES****/

/*get counts*/
frame create stdizemortallbycause
frame change stdizemortallbycause

use obj1_deaths_2010_2019

collapse (count) total_mortality_count=_temp, by(year agegroup_5yr underlying_cause_cat)   

drop if underlying_cause_cat==. | agegroup_5yr==.

/*generate new rows for causes of death with 0 counts for a given age group - as these are missing from the dataset after the above collapse (due to 0 count)*/
fillin agegroup_5yr year underlying_cause_cat

replace total_mortality_count = 0 if missing(total_mortality_count) /*replace missing counts for the var total_mortality_count for the newly created rows with 0*/

/*merge in population data*/
frlink m:1 year agegroup_5yr, frame(population5yrbands)

frget population, from(population5yrbands)

drop population5yrbands

/*merge in WHO standard population weights*/
frlink m:1 agegroup_5yr, frame(weightsall)

frget weight_whostd_all, from(weightsall)

drop weightsall

/*calculate standardized rates*/
generate crude = total_mortality_count / population

generate product = crude*weight_whostd_all
gen var_component = (weight_whostd_all^2) * (total_mortality_count / (population^2))

egen _total_mortality_rate_adj = total(product), by(year underlying_cause_cat)
egen _se_adj = total(var_component), by(year underlying_cause_cat)

generate total_mortality_rate_adj = _total_mortality_rate_adj*100000
gen se_adj=sqrt(_se_adj) * 100000

contract year underlying_cause_cat total_mortality_rate_adj se_adj

format total_mortality_rate_adj %9.0gc

drop _freq


/*reshape and rename variables*/

reshape wide total_mortality_rate_adj se_adj, i(year) j(underlying_cause_cat)

local codcat "matperi infect cancer diabendo cvd resp digest gu psyc external other"
local j = 1
forvalues j = 1/11 {
	local sffx : word `j' of `codcat' 
	rename (total_mortality_rate_adj`j' se_adj`j') (total_mortality_rateadj_`sffx' se_adj_`sffx')
    local ++j
}


*Obtain cause-specific age-standardized rates of ypll 70*/

/*get counts*/
frame create stdizeypllallcause
frame change stdizeypllallcause

use obj1_deaths_2010_2019

collapse ///
    (sum) total_ypll_count = ypll_70_burden ///
	(sum) total_ypll_squared = ypll_70_sq_burden, ///
    by(year agegroup_5yr underlying_cause_cat)

drop if underlying_cause_cat==. | agegroup_5yr==.

/*generate new rows for causes of death with 0 counts for a given age group - as these are missing from the dataset after the above collapse (due to 0 count)*/
fillin agegroup_5yr year underlying_cause_cat

/*replace missing counts for the newly created rows with 0*/
replace total_ypll_count = 0 if missing(total_ypll_count)

replace total_ypll_squared = 0 if missing(total_ypll_squared)


/*merge in population data*/
frlink m:1 year agegroup_5yr, frame(population5yrbands)
frget population, from(population5yrbands)
drop population5yrbands

/*merge in WHO standard population weights (for subset of pop <70 years)*/
frlink m:1 agegroup_5yr, frame(weights70)
frget weight_whostd_70, from(weights70)
drop weights70

tab1 agegroup_5yr, nolabel 
drop if agegroup_5yr >= 15 

/*calculate standardized rates*/
generate crude = total_ypll_count / population
generate product = crude*weight_whostd_70

gen double var_component = (weight_whostd_70^2 * total_ypll_squared) / (population^2)

collapse ///
    (sum) _total_ypll_rate_adj = product ///
    (sum) _se_ypll_adj = var_component, ///
  by(year underlying_cause_cat)

generate total_ypll_rate_adj = _total_ypll_rate_adj * 100000
gen se_ypll_adj=(sqrt(_se_ypll_adj)) * 100000

format total_ypll_rate_adj se_ypll_adj %9.0gc


/*reshape and rename variables*/

drop _total_ypll_rate_adj _se_ypll_adj

reshape wide total_ypll_rate_adj se_ypll_adj, i(year) j(underlying_cause_cat)

local codcat "matperi infect cancer diabendo cvd resp digest gu psyc external other"
local j = 1
forvalues j = 1/11 {
	local sffx : word `j' of `codcat' 
	rename (total_ypll_rate_adj`j' se_ypll_adj`j') (total_ypll_rateadj_`sffx' se_ypll_adj_`sffx')
    local ++j
}

/*link frames with the adjusted rates and output as Excel file*/

*Total mortality - all cause

frlink 1:1 year, frame(stdizemortall)

frget total_mortality_rate_adj_all se_adj_all, from(stdizemortall)

drop stdizemortall

*Premature mortality - all cause

frlink 1:1 year, frame(stdizeypllall)

frget total_ypll_rate_adj_all se_ypll_adj_all, from(stdizeypllall)

drop stdizeypllall

*Total mortality - by cause

frlink 1:1 year, frame(stdizemortallbycause)

frget total_mortality_rateadj_matperi total_mortality_rateadj_infect total_mortality_rateadj_cancer total_mortality_rateadj_diabendo total_mortality_rateadj_cvd total_mortality_rateadj_resp total_mortality_rateadj_digest total_mortality_rateadj_gu total_mortality_rateadj_psyc total_mortality_rateadj_external total_mortality_rateadj_other se_adj_matperi se_adj_infect se_adj_cancer se_adj_diabendo se_adj_cvd se_adj_resp se_adj_digest se_adj_gu se_adj_psyc se_adj_external se_adj_other, from(stdizemortallbycause)

drop stdizemortallbycause

label variable total_mortality_rate_adj_all "Total Mortality Rate (Adjusted to WHO Std) - Overall"
label variable total_mortality_rateadj_matperi "Total Mortality Rate (Adjusted to WHO Std) - Maternal & Perinatal Conditions"
label variable total_mortality_rateadj_infect "Total Mortality Rate (Adjusted to WHO Std) - Infectious Diseases"
label variable total_mortality_rateadj_cancer "Total Mortality Rate (Adjusted to WHO Std) - Cancer"
label variable total_mortality_rateadj_diabendo "Total Mortality Rate (Adjusted to WHO Std) - Diabetes & Endocrine Disorders"
label variable total_mortality_rateadj_cvd "Total Mortality Rate (Adjusted to WHO Std) - Circulatory Diseases"
label variable total_mortality_rateadj_resp "Total Mortality Rate (Adjusted to WHO Std) - Respiratory Diseases"
label variable total_mortality_rateadj_digest "Total Mortality Rate (Adjusted to WHO Std) - Digestive Diseases"
label variable total_mortality_rateadj_gu "Total Mortality Rate (Adjusted to WHO Std) - Genitourinary Diseases"
label variable total_mortality_rateadj_psyc "Total Mortality Rate (Adjusted to WHO Std) - Neuropsychiatric Conditions"
label variable total_mortality_rateadj_external "Total Mortality Rate (Adjusted to WHO Std) - External Causes"
label variable total_mortality_rateadj_other "Total Mortality Rate (Adjusted to WHO Std) - Other"

label variable total_ypll_rate_adj_all "Total YPLL Rate (Adjusted to WHO Std) - Overall"
label variable total_ypll_rateadj_matperi "Total YPLL Rate (Adjusted to WHO Std) - Maternal & Perinatal Conditions"
label variable total_ypll_rateadj_infect "Total YPLL Rate (Adjusted to WHO Std) - Infectious Diseases"
label variable total_ypll_rateadj_cancer "Total YPLL Rate (Adjusted to WHO Std) - Cancer"
label variable total_ypll_rateadj_diabendo "Total YPLL Rate (Adjusted to WHO Std) - Diabetes & Endocrine Disorders"
label variable total_ypll_rateadj_cvd "Total YPLL Rate (Adjusted to WHO Std) - Circulatory Diseases"
label variable total_ypll_rateadj_resp "Total YPLL Rate (Adjusted to WHO Std) - Respiratory Diseases"
label variable total_ypll_rateadj_digest "Total YPLL Rate (Adjusted to WHO Std) - Digestive Diseases"
label variable total_ypll_rateadj_gu "Total YPLL Rate (Adjusted to WHO Std) - Genitourinary Diseases"
label variable total_ypll_rateadj_psyc "Total YPLL Rate (Adjusted to WHO Std) - Neuropsychiatric Conditions"
label variable total_ypll_rateadj_external "Total YPLL Rate (Adjusted to WHO Std) - External Causes"
label variable total_ypll_rateadj_other "Total YPLL Rate (Adjusted to WHO Std) - Other"

label variable se_adj_all "Standard Error for Adjusted Mortality Rate - Overall"
label variable se_adj_matperi "Standard Error for Adjusted Mortality Rate - Maternal & Perinatal Conditions"
label variable se_adj_infect "Standard Error for Adjusted Mortality Rate - Infectious Diseases"
label variable se_adj_cancer "Standard Error for Adjusted Mortality Rate - Cancer"
label variable se_adj_diabendo "Standard Error for Adjusted Mortality Rate - Diabetes & Endocrine Disorders"
label variable se_adj_cvd "Standard Error for Adjusted Mortality Rate - Circulatory Diseases"
label variable se_adj_resp "Standard Error for Adjusted Mortality Rate - Respiratory Diseases"
label variable se_adj_digest "Standard Error for Adjusted Mortality Rate - Digestive Diseases"
label variable se_adj_gu "Standard Error for Adjusted Mortality Rate - Genitourinary Diseases"
label variable se_adj_psyc "Standard Error for Adjusted Mortality Rate - Neuropsychiatric Conditions"
label variable se_adj_external "Standard Error for Adjusted Mortality Rate - External Causes"
label variable se_adj_other "Standard Error for Adjusted Mortality Rate - Other"

label variable se_ypll_adj_all "Standard Error for Adjusted YPLL Rate - Overall"
label variable se_ypll_adj_matperi "Standard Error for Adjusted YPLL Rate - Maternal & Perinatal Conditions"
label variable se_ypll_adj_infect "Standard Error for Adjusted YPLL Rate - Infectious Diseases"
label variable se_ypll_adj_cancer "Standard Error for Adjusted YPLL Rate - Cancer"
label variable se_ypll_adj_diabendo "Standard Error for Adjusted YPLL Rate - Diabetes & Endocrine Disorders"
label variable se_ypll_adj_cvd "Standard Error for Adjusted YPLL Rate - Circulatory Diseases"
label variable se_ypll_adj_resp "Standard Error for Adjusted YPLL Rate - Respiratory Diseases"
label variable se_ypll_adj_digest "Standard Error for Adjusted YPLL Rate - Digestive Diseases"
label variable se_ypll_adj_gu "Standard Error for Adjusted YPLL Rate - Genitourinary Diseases"
label variable se_ypll_adj_psyc "Standard Error for Adjusted YPLL Rate - Neuropsychiatric Conditions"
label variable se_ypll_adj_external "Standard Error for Adjusted YPLL Rate - External Causes"
label variable se_ypll_adj_other "Standard Error for Adjusted YPLL Rate - Other"

order year total_mortality_rate_adj_all total_mortality_rateadj_matperi total_mortality_rateadj_infect total_mortality_rateadj_cancer  total_mortality_rateadj_diabendo total_mortality_rateadj_cvd total_mortality_rateadj_resp  total_mortality_rateadj_digest total_mortality_rateadj_gu total_mortality_rateadj_psyc total_mortality_rateadj_external total_mortality_rateadj_other total_ypll_rate_adj_all total_ypll_rateadj_matperi total_ypll_rateadj_infect total_ypll_rateadj_cancer total_ypll_rateadj_diabendo total_ypll_rateadj_cvd total_ypll_rateadj_resp total_ypll_rateadj_digest total_ypll_rateadj_gu total_ypll_rateadj_psyc total_ypll_rateadj_external total_ypll_rateadj_other se_adj_all se_adj_matperi se_adj_infect se_adj_cancer se_adj_diabendo se_adj_cvd se_adj_resp se_adj_digest se_adj_gu se_adj_psyc se_adj_external se_adj_other se_ypll_adj_all se_ypll_adj_matperi se_ypll_adj_infect se_ypll_adj_cancer se_ypll_adj_diabendo se_ypll_adj_cvd se_ypll_adj_resp se_ypll_adj_digest se_ypll_adj_gu se_ypll_adj_psyc se_ypll_adj_external se_ypll_adj_other

save obj1_overall70_stdize, replace





/************************
 SECTION 2 - BY SEX
************************/

/****SECTION 2.1 - BY SEX, ALL CAUSE MORTALITY - CRUDE RATES****/

frame create sexoverall
frame change sexoverall

use obj1_deaths_2010_2019

collapse /// 
	(count) total_mortality_count = _temp /// 
	(mean) premature_percentage = premature70 /// 
	(sum) premature_count = premature70 ///
	(sum) total_ypll_count = ypll_70 ///
	(sum) total_ypll_squared = ypll_70_sq_burden ///
	(median) median_ypll = ypll_70 ///
	(mean) mean_ypll = ypll_70, ///
by(year sex)

replace premature_percentage = premature_percentage * 100

drop if !inlist(sex, 1, 2)

format total_mortality_count total_ypll_count premature_count %9.0gc
format premature_percentage %8.2f
format median_ypll mean_ypll %4.2f


/*Add in population data*/

*overall population

frame create populationdatasex
frame change populationdatasex

use "${clean_data_dir}/population_sex.dta" 

frame change sexoverall

frlink 1:1 year sex, frame(populationdatasex)

frget population, from(populationdatasex)

*population under age 70

frame create populationdata70sex
frame change populationdata70sex

use "${clean_data_dir}/population_agesex.dta" 

tab1 agegroup

drop if inlist(agegroup,7,8)

tab1 agegroup

collapse (sum) population70=population, by(year sex) 

frame change sexoverall

frlink 1:1 year sex, frame(populationdata70sex)

frget population70, from(populationdata70sex)

/*Calculate rates*/

generate total_mortality_rate = (total_mortality_count/population)*100000
generate total_ypll_rate = (total_ypll_count/population70)*100000
generate se_ypll_rate = (sqrt(total_ypll_squared)/population70)*100000

drop total_ypll_squared

format total_mortality_rate total_ypll_rate se_ypll_rate %9.0gc

order year sex population total_mortality_count total_mortality_rate premature_count premature_percentage total_ypll_count total_ypll_rate se_ypll_rate median_ypll mean_ypll

rename (total_mortality_count total_mortality_rate premature_count premature_percentage total_ypll_count total_ypll_rate se_ypll_rate median_ypll mean_ypll) (total_mortality_count_all total_mortality_rate_all premature_count_all premature_percentage_all total_ypll_count_all total_ypll_rate_all se_ypll_rate_all median_ypll_all mean_ypll_all)

drop populationdatasex populationdata70sex


/*indicator 5 - unconditional probability of death between ages 30 and 70 from NCDs*/
frame create ncdprobsex
frame change ncdprobsex

use obj1_deaths_2010_2019

collapse (sum) total_mortality_count_ncd=underlying_cause_sdgncd, by(year agegroup_5yr sex)

tab1 agegroup, missing
tab1 agegroup, nolabel 

drop if agegroup==. | agegroup<7 | agegroup>14 /*probability to be calculated between ages 30 and 70*/

tab1 agegroup, missing

tab1 sex, missing /* no unknown values */

/*merge in population data*/
frame create population5yrbandssex
frame change population5yrbandssex

use "${clean_data_dir}/population_agesex_5yrbands.dta" 

frame change ncdprobsex

frlink 1:1 year agegroup_5yr sex, frame(population5yrbandssex)

frget population, from(population5yrbandssex)

drop population5yrbandssex

/*calculate unconditional probability of dying*/
generate mortrate = total_mortality_count_ncd / population /*step 1: calculate age-specific mortality rates for each five-year age group*/

generate mortprob = mortrate*5 / (1+mortrate*2.5) /*step 2: translate the mortality rate into the probability of death in each 5-yer age range*/

/*step 3: multiply the probabilities across age groups to get the total probability*/
keep agegroup_5yr sex year mortprob total_mortality_count_ncd

reshape wide total_mortality_count_ncd mortprob, i(year sex) j(agegroup_5yr)

generate probability_death_ncd = 1-((1-mortprob7)*(1-mortprob8)*(1-mortprob9)*(1-mortprob10)*(1-mortprob11)*(1-mortprob12)*(1-mortprob13)*(1-mortprob14))

/*obtain the count to correspond with the percentage*/
generate mortality_count_ncd = total_mortality_count_ncd7 + total_mortality_count_ncd8 + total_mortality_count_ncd9 + total_mortality_count_ncd10 + total_mortality_count_ncd11 + total_mortality_count_ncd12 + total_mortality_count_ncd13 + total_mortality_count_ncd14

format probability_death_ncd %04.3f
format mortality_count_ncd %9.0gc

/*link back to original frame*/

frame change sexoverall

frlink 1:1 year sex, frame(ncdprobsex)

frget mortality_count_ncd probability_death_ncd, from(ncdprobsex)

drop ncdprobsex


/****SECTION 2.2 - BY SEX, ALL CAUSE MORTALITY - AGE-STANDARDIZED RATES****/

/*Total Mortality*/

/*get counts*/
frame create stdizemortsex
frame change stdizemortsex

use obj1_deaths_2010_2019

collapse (count) total_mortality_count=_temp, by(year agegroup_5yr sex)   

drop if agegroup==. | sex==9 /*no filling in needed -- all combos of age and sex and year represented*/

/*merge in population data*/
frlink 1:1 year agegroup_5yr sex, frame(population5yrbandssex)

frget population, from(population5yrbandssex)

drop population5yrbandssex

/*merge in WHO standard population weights*/
frlink m:1 agegroup_5yr, frame(weightsall)

frget weight_whostd_all, from(weightsall)

drop weightsall


/*calculate standardized rates*/
generate crude = total_mortality_count / population

generate product = crude*weight_whostd_all
gen var_component = (weight_whostd_all^2) * (total_mortality_count / (population^2))

egen _total_mortality_rate_adj_all = total(product), by(year sex)
egen _se_adj_all = total(var_component), by(year sex)

generate total_mortality_rate_adj_all = _total_mortality_rate_adj_all*100000
gen se_adj_all=sqrt(_se_adj_all) * 100000

contract year sex total_mortality_rate_adj_all se_adj_all

format total_mortality_rate_adj_all %9.0gc

drop _freq


/*Premature Mortality*/

/*get counts*/
frame create stdizeypllsex
frame change stdizeypllsex

use obj1_deaths_2010_2019

collapse ///
    (sum) total_ypll_count = ypll_70_burden ///
	(sum) total_ypll_squared = ypll_70_sq_burden, ///
    by(year agegroup_5yr sex)

drop if agegroup_5yr==. | sex==9 /*no filling in needed -- all combos of age and sex and year represented*/

/*merge in population data*/
frlink 1:1 year agegroup_5yr sex, frame(population5yrbandssex)

frget population, from(population5yrbandssex)

drop population5yrbandssex

/*merge in WHO standard population weights (for subset of pop <70 years)*/
frlink m:1 agegroup_5yr, frame(weights70)
frget weight_whostd_70, from(weights70)
drop weights70

tab1 agegroup_5yr, nolabel 
drop if agegroup_5yr >= 15 

/*calculate standardized rates*/
generate crude = total_ypll_count / population
generate product = crude*weight_whostd_70

gen double var_component = (weight_whostd_70^2 * total_ypll_squared) / (population^2)

collapse ///
    (sum) _total_ypll_rate_adj_all = product ///
    (sum) _se_ypll_adj_all = var_component, ///
    by(year sex)

generate total_ypll_rate_adj_all = _total_ypll_rate_adj_all * 100000
gen se_ypll_adj_all = (sqrt(_se_ypll_adj_all)) * 100000

format total_ypll_rate_adj_all se_ypll_adj_all %9.0gc


/****SECTION 2.3 - BY SEX STRATIFIED BY CAUSE - CRUDE RATES****/

frame create sexbymjrcause
frame change sexbymjrcause

use obj1_deaths_2010_2019

collapse ///
	(count) total_mortality_count = _temp ///
	(mean) premature_percentage = premature70 ///
	(sum) premature_count = premature70 ///
	(sum) total_ypll_count = ypll_70 ///
	(sum) total_ypll_squared = ypll_70_sq_burden ///
	(median) median_ypll = ypll_70 ///
	(mean) mean_ypll = ypll_70, ///
by(year sex underlying_cause_cat) 

replace premature_percentage = premature_percentage * 100

drop if underlying_cause_cat==. | sex==9

format total_mortality_count total_ypll_count premature_count %9.0gc
format premature_percentage %8.2f
format median_ypll mean_ypll %4.2f

reshape wide total_mortality_count total_ypll_count total_ypll_squared premature_count premature_percentage median_ypll mean_ypll, i(year sex) j(underlying_cause_cat)

/*Add in population data and calculate rates*/

*overall population

frlink 1:1 year sex, frame(populationdatasex)

frget population, from(populationdatasex)

*population under age 70

frlink 1:1 year sex, frame(populationdata70sex)

frget population70, from(populationdata70sex)


*calculate rates of mortality and premature mortality, rename and label variables

local codcat "matperi infect cancer diabendo cvd resp digest gu psyc external other"
local j = 1
forvalues j = 1/11 {
	local sffx : word `j' of `codcat' 
	rename (total_mortality_count`j' premature_percentage`j' premature_count`j' total_ypll_count`j' total_ypll_squared`j' median_ypll`j' mean_ypll`j') (total_mortality_count_`sffx' premature_percentage_`sffx' premature_count_`sffx' total_ypll_count_`sffx' total_ypll_squared_`sffx' median_ypll_`sffx' mean_ypll_`sffx')
    local ++j
}

local codcat2 "matperi infect cancer diabendo cvd resp digest gu psyc external other"
local q = 1
foreach x of local codcat2 {
    generate total_mortality_rate_`x' = (total_mortality_count_`x'/population)*100000
	generate total_ypll_rate_`x' = (total_ypll_count_`x'/population70)*100000
	generate se_ypll_rate_`x' = (sqrt(total_ypll_squared_`x')/population70)*100000
	format total_mortality_rate_`x' total_ypll_rate_`x' se_ypll_rate_`x' %9.0gc
    local ++q
}

drop total_ypll_squared_*

/*link frames with crude rates, label variables*/

frlink 1:1 year sex, frame(sexoverall)

frget total_mortality_count_all total_mortality_rate_all premature_count_all premature_percentage_all total_ypll_count_all total_ypll_rate_all se_ypll_rate_all median_ypll_all mean_ypll_all mortality_count_ncd probability_death_ncd, from(sexoverall)

label variable population70 "Jamaica Mid-Year Population Estimates (Pop <70 Years of Age)"

label variable total_mortality_count_all "Total Mortality Count - Overall"
label variable total_mortality_rate_all "Total Mortality Rate - Overall"
label variable premature_count_all "No. Premature Deaths - Overall"
label variable premature_percentage_all "Pctg Premature Deaths - Overall"
label variable total_ypll_count_all "Total YPLL Count - Overall"
label variable total_ypll_rate_all "Total YPLL Rate - Overall"
label variable median_ypll_all "Median YPLL/Person - Overall"
label variable mean_ypll_all "Mean YPLL/Person - Overall"
label variable mortality_count_ncd "Number of Premature NCD Deaths (WHO Defn)"
label variable probability_death_ncd "Probability of Premature NCD Death (WHO Defn)"

label variable total_mortality_count_matperi "Total Mortality Count - Maternal & Perinatal Conditions"
label variable total_mortality_rate_matperi "Total Mortality Rate - Maternal & Perinatal Conditions"
label variable premature_count_matperi "No. Premature Deaths - Maternal & Perinatal Conditions"
label variable premature_percentage_matperi "Pctg Premature Deaths - Maternal & Perinatal Conditions"
label variable total_ypll_count_matperi "Total YPLL Count - Maternal & Perinatal Conditions"
label variable total_ypll_rate_matperi "Total YPLL Rate - Maternal & Perinatal Conditions"
label variable median_ypll_matperi "Median YPLL/Person - Maternal & Perinatal Conditions"
label variable mean_ypll_matperi "Mean YPLL/Person - Maternal & Perinatal Conditions"

label variable total_mortality_count_infect "Total Mortality Count - Infectious Diseases"
label variable total_mortality_rate_infect "Total Mortality Rate - Infectious Diseases"
label variable premature_count_infect "No. Premature Deaths - Infectious Diseases"
label variable premature_percentage_infect "Pctg Premature Deaths - Infectious Diseases"
label variable total_ypll_count_infect "Total YPLL Count - Infectious Diseases"
label variable total_ypll_rate_infect "Total YPLL Rate - Infectious Diseases"
label variable median_ypll_infect "Median YPLL/Person - Infectious Diseases"
label variable mean_ypll_infect "Mean YPLL/Person - Infectious Diseases"

label variable total_mortality_count_cancer "Total Mortality Count - Cancer"
label variable total_mortality_rate_cancer "Total Mortality Rate - Cancer"
label variable premature_count_cancer "No. Premature Deaths - Cancer"
label variable premature_percentage_cancer "Pctg Premature Deaths - Cancer"
label variable total_ypll_count_cancer "Total YPLL Count - Cancer"
label variable total_ypll_rate_cancer "Total YPLL Rate - Cancer"
label variable median_ypll_cancer "Median YPLL/Person - Cancer"
label variable mean_ypll_cancer "Mean YPLL/Person - Cancer"

label variable total_mortality_count_diabendo "Total Mortality Count - Diabetes & Endocrine Disorders"
label variable total_mortality_rate_diabendo "Total Mortality Rate - Diabetes & Endocrine Disorders"
label variable premature_count_diabendo "No. Premature Deaths - Diabetes & Endocrine Disorders"
label variable premature_percentage_diabendo "Pctg Premature Deaths - Diabetes & Endocrine Disorders"
label variable total_ypll_count_diabendo "Total YPLL Count - Diabetes & Endocrine Disorders"
label variable total_ypll_rate_diabendo "Total YPLL Rate - Diabetes & Endocrine Disorders"
label variable median_ypll_diabendo "Median YPLL/Person - Diabetes & Endocrine Disorders"
label variable mean_ypll_diabendo "Mean YPLL/Person - Diabetes & Endocrine Disorders"

label variable total_mortality_count_cvd "Total Mortality Count - Circulatory Diseases"
label variable total_mortality_rate_cvd "Total Mortality Rate - Circulatory Diseases"
label variable premature_count_cvd "No. Premature Deaths - Circulatory Diseases"
label variable premature_percentage_cvd "Pctg Premature Deaths - Circulatory Diseases"
label variable total_ypll_count_cvd "Total YPLL Count - Circulatory Diseases"
label variable total_ypll_rate_cvd "Total YPLL Rate - Circulatory Diseases"
label variable median_ypll_cvd "Median YPLL/Person - Circulatory Diseases"
label variable mean_ypll_cvd "Mean YPLL/Person - Circulatory Diseases"

label variable total_mortality_count_resp "Total Mortality Count - Respiratory Diseases"
label variable total_mortality_rate_resp "Total Mortality Rate - Respiratory Diseases"
label variable premature_count_resp "No. Premature Deaths - Respiratory Diseases"
label variable premature_percentage_resp "Pctg Premature Deaths - Respiratory Diseases"
label variable total_ypll_count_resp "Total YPLL Count - Respiratory Diseases"
label variable total_ypll_rate_resp "Total YPLL Rate - Respiratory Diseases"
label variable median_ypll_resp "Median YPLL/Person - Respiratory Diseases"
label variable mean_ypll_resp "Mean YPLL/Person - Respiratory Diseases"

label variable total_mortality_count_digest "Total Mortality Count - Digestive Diseases"
label variable total_mortality_rate_digest "Total Mortality Rate - Digestive Diseases"
label variable premature_count_digest "No. Premature Deaths - Digestive Diseases"
label variable premature_percentage_digest "Pctg Premature Deaths - Digestive Diseases"
label variable total_ypll_count_digest "Total YPLL Count - Digestive Diseases"
label variable total_ypll_rate_digest "Total YPLL Rate - Digestive Diseases"
label variable median_ypll_digest "Median YPLL/Person - Digestive Diseases"
label variable mean_ypll_digest "Mean YPLL/Person - Digestive Diseases"

label variable total_mortality_count_gu "Total Mortality Count - Genitourinary Diseases"
label variable total_mortality_rate_gu "Total Mortality Rate - Genitourinary Diseases"
label variable premature_count_gu "No. Premature Deaths - Genitourinary Diseases"
label variable premature_percentage_gu "Pctg Premature Deaths - Genitourinary Diseases"
label variable total_ypll_count_gu "Total YPLL Count - Genitourinary Diseases"
label variable total_ypll_rate_gu "Total YPLL Rate - Genitourinary Diseases"
label variable median_ypll_gu "Median YPLL/Person - Genitourinary Diseases"
label variable mean_ypll_gu "Mean YPLL/Person - Genitourinary Diseases"

label variable total_mortality_count_psyc "Total Mortality Count - Neuropsychiatric Conditions"
label variable total_mortality_rate_psyc "Total Mortality Rate - Neuropsychiatric Conditions"
label variable premature_count_psyc "No. Premature Deaths - Neuropsychiatric Conditions"
label variable premature_percentage_psyc "Pctg Premature Deaths - Neuropsychiatric Conditions"
label variable total_ypll_count_psyc "Total YPLL Count - Neuropsychiatric Conditions"
label variable total_ypll_rate_psyc "Total YPLL Rate - Neuropsychiatric Conditions"
label variable median_ypll_psyc "Median YPLL/Person - Neuropsychiatric Conditions"
label variable mean_ypll_psyc "Mean YPLL/Person - Neuropsychiatric Conditions"

label variable total_mortality_count_external "Total Mortality Count - External Causes"
label variable total_mortality_rate_external "Total Mortality Rate - External Causes"
label variable premature_count_external "No. Premature Deaths - External Causes"
label variable premature_percentage_external "Pctg Premature Deaths - External Causes"
label variable total_ypll_count_external "Total YPLL Count - External Causes"
label variable total_ypll_rate_external "Total YPLL Rate - External Causes"
label variable median_ypll_external "Median YPLL/Person - External Causes"
label variable mean_ypll_external "Mean YPLL/Person - External Causes"

label variable total_mortality_count_other "Total Mortality Count - Other"
label variable total_mortality_rate_other "Total Mortality Rate - Other"
label variable premature_count_other "No. Premature Deaths - Other"
label variable premature_percentage_other "Pctg Premature Deaths - Other"
label variable total_ypll_count_other "Total YPLL Count - Other"
label variable total_ypll_rate_other "Total YPLL Rate - Other"
label variable median_ypll_other "Median YPLL/Person - Other"
label variable mean_ypll_other "Mean YPLL/Person - Other"

label variable se_ypll_rate_all "Standard Error for Crude YPLL Rate - Overall"
label variable se_ypll_rate_matperi "Standard Error for Crude YPLL Rate - Maternal & Perinatal Conditions"
label variable se_ypll_rate_infect "Standard Error for Crude YPLL Rate - Infectious Diseases"
label variable se_ypll_rate_cancer "Standard Error for Crude YPLL Rate - Cancer"
label variable se_ypll_rate_diabendo "Standard Error for Crude YPLL Rate - Diabetes & Endocrine Disorders"
label variable se_ypll_rate_cvd "Standard Error for Crude YPLL Rate - Circulatory Diseases"
label variable se_ypll_rate_resp "Standard Error for Crude YPLL Rate - Respiratory Diseases"
label variable se_ypll_rate_digest "Standard Error for Crude YPLL Rate - Digestive Diseases"
label variable se_ypll_rate_gu "Standard Error for Crude YPLL Rate - Genitourinary Diseases"
label variable se_ypll_rate_psyc "Standard Error for Crude YPLL Rate - Neuropsychiatric Conditions"
label variable se_ypll_rate_external "Standard Error for Crude YPLL Rate - External Causes"
label variable se_ypll_rate_other "Standard Error for Crude YPLL Rate - Other"


order year sex population population70 total_mortality_count_all total_mortality_rate_all total_mortality_count_matperi total_mortality_rate_matperi total_mortality_count_infect total_mortality_rate_infect total_mortality_count_cancer total_mortality_rate_cancer total_mortality_count_diabendo total_mortality_rate_diabendo total_mortality_count_cvd total_mortality_rate_cvd total_mortality_count_resp total_mortality_rate_resp total_mortality_count_digest total_mortality_rate_digest total_mortality_count_gu total_mortality_rate_gu total_mortality_count_psyc total_mortality_rate_psyc total_mortality_count_external total_mortality_rate_external total_mortality_count_other total_mortality_rate_other premature_count_all premature_percentage_all premature_count_matperi premature_percentage_matperi premature_count_infect premature_percentage_infect premature_count_cancer premature_percentage_cancer premature_count_diabendo premature_percentage_diabendo premature_count_cvd premature_percentage_cvd premature_count_resp premature_percentage_resp premature_count_digest premature_percentage_digest premature_count_gu premature_percentage_gu premature_count_psyc premature_percentage_psyc premature_count_external premature_percentage_external premature_count_other premature_percentage_other total_ypll_count_all total_ypll_rate_all total_ypll_count_matperi total_ypll_rate_matperi total_ypll_count_infect total_ypll_rate_infect total_ypll_count_cancer total_ypll_rate_cancer total_ypll_count_diabendo total_ypll_rate_diabendo total_ypll_count_cvd total_ypll_rate_cvd total_ypll_count_resp total_ypll_rate_resp total_ypll_count_digest total_ypll_rate_digest total_ypll_count_gu total_ypll_rate_gu total_ypll_count_psyc total_ypll_rate_psyc total_ypll_count_external total_ypll_rate_external total_ypll_count_other total_ypll_rate_other median_ypll_all median_ypll_matperi median_ypll_infect median_ypll_cancer median_ypll_diabendo median_ypll_cvd median_ypll_resp median_ypll_digest median_ypll_gu median_ypll_psyc median_ypll_external median_ypll_other mean_ypll_all mean_ypll_matperi mean_ypll_infect mean_ypll_cancer mean_ypll_diabendo mean_ypll_cvd mean_ypll_resp mean_ypll_digest mean_ypll_gu mean_ypll_psyc mean_ypll_external mean_ypll_other mortality_count_ncd probability_death_ncd se_ypll_rate_all se_ypll_rate_matperi se_ypll_rate_infect se_ypll_rate_cancer se_ypll_rate_diabendo se_ypll_rate_cvd se_ypll_rate_resp se_ypll_rate_digest se_ypll_rate_gu se_ypll_rate_psyc se_ypll_rate_external se_ypll_rate_other

sort year sex

drop populationdatasex populationdata70sex sexoverall 

save obj1_sex70, replace


/****SECTION 2.4 - BY SEX STRATIFIED BY CAUSE - AGE-STANDARDIZED RATES****/

/*get counts*/
frame create stdizemortsexbycause
frame change stdizemortsexbycause

use obj1_deaths_2010_2019

collapse (count) total_mortality_count=_temp, by(year agegroup_5yr sex underlying_cause_cat)   

drop if underlying_cause_cat==. | agegroup_5yr==. | sex==9

/*generate new rows for causes of death with 0 counts for a given sex and age group combo - as these are missing from the dataset after the above collapse (due to 0 count)*/
fillin agegroup_5yr sex year underlying_cause_cat

replace total_mortality_count = 0 if missing(total_mortality_count) /*replace missing counts for the var total_mortality_count for the newly created rows with 0*/

/*merge in population data*/
frlink m:1 year agegroup_5yr sex, frame(population5yrbandssex)

frget population, from(population5yrbandssex)

drop population5yrbandssex

/*merge in WHO standard population weights*/
frlink m:1 agegroup_5yr, frame(weightsall)

frget weight_whostd_all, from(weightsall)

drop weightsall

/*calculate standardized rates*/
generate crude = total_mortality_count / population

generate product = crude*weight_whostd_all
gen var_component = (weight_whostd_all^2) * (total_mortality_count / (population^2))

egen _total_mortality_rate_adj = total(product), by(year sex underlying_cause_cat)
egen _se_adj 				   = total(var_component), by(year sex underlying_cause_cat)

generate total_mortality_rate_adj = _total_mortality_rate_adj*100000
gen se_adj=sqrt(_se_adj) * 100000

contract year sex underlying_cause_cat total_mortality_rate_adj se_adj

format total_mortality_rate_adj %9.0gc

drop _freq

/*reshape and rename variables*/

reshape wide total_mortality_rate_adj se_adj, i(year sex) j(underlying_cause_cat)

local codcat "matperi infect cancer diabendo cvd resp digest gu psyc external other"
local j = 1
forvalues j = 1/11 {
	local sffx : word `j' of `codcat' 
	rename (total_mortality_rate_adj`j' se_adj`j') (total_mortality_rateadj_`sffx' se_adj_`sffx')
    local ++j
}

/*Obtain cause-specific age-standardized rates of ypll 70*/

/*get counts*/
frame create stdizeypllsexcause
frame change stdizeypllsexcause

use obj1_deaths_2010_2019

collapse ///
    (sum) total_ypll_count = ypll_70_burden ///
	(sum) total_ypll_squared = ypll_70_sq_burden, ///
    by(year agegroup_5yr sex underlying_cause_cat)

drop if underlying_cause_cat==. | agegroup_5yr==. | sex==9

/*generate new rows for causes of death with 0 counts for a given sex and age group combo - as these are missing from the dataset after the above collapse (due to 0 count)*/
fillin agegroup_5yr sex year underlying_cause_cat

 /*replace missing counts for the newly created rows with 0*/
replace total_ypll_count = 0 if missing(total_ypll_count)

replace total_ypll_squared = 0 if missing(total_ypll_squared)


/*merge in population data*/
frlink m:1 year agegroup_5yr sex, frame(population5yrbandssex)
frget population, from(population5yrbandssex)
drop population5yrbandssex

/*merge in WHO standard population weights (for subset of pop <70 years)*/
frlink m:1 agegroup_5yr, frame(weights70)
frget weight_whostd_70, from(weights70)
drop weights70

tab1 agegroup_5yr, nolabel 
drop if agegroup_5yr >= 15 

/*calculate standardized rates*/
generate crude = total_ypll_count / population
generate product = crude*weight_whostd_70

gen double var_component = (weight_whostd_70^2 * total_ypll_squared) / (population^2)

collapse ///
    (sum) _total_ypll_rate_adj = product ///
    (sum) _se_ypll_adj = var_component, ///
    by(year sex underlying_cause_cat)

generate total_ypll_rate_adj = _total_ypll_rate_adj * 100000
gen se_ypll_adj=(sqrt(_se_ypll_adj)) * 100000

format total_ypll_rate_adj se_ypll_adj %9.0gc

/*reshape and rename variables*/

drop _total_ypll_rate_adj _se_ypll_adj

reshape wide total_ypll_rate_adj se_ypll_adj, i(year sex) j(underlying_cause_cat)

local codcat "matperi infect cancer diabendo cvd resp digest gu psyc external other"
local j = 1
forvalues j = 1/11 {
	local sffx : word `j' of `codcat' 
	rename (total_ypll_rate_adj`j' se_ypll_adj`j') (total_ypll_rateadj_`sffx' se_ypll_adj_`sffx')
    local ++j
}

/*link frames with the adjusted rates*/

*Total mortality - all cause

frlink 1:1 year sex, frame(stdizemortsex)

frget total_mortality_rate_adj_all se_adj_all, from(stdizemortsex)

drop stdizemortsex

*Premature mortality - all cause

frlink 1:1 year sex, frame(stdizeypllsex)

frget total_ypll_rate_adj_all se_ypll_adj_all, from(stdizeypllsex)

drop stdizeypllsex

*Total mortality - by cause

frlink 1:1 year sex, frame(stdizemortsexbycause)

frget total_mortality_rateadj_matperi total_mortality_rateadj_infect total_mortality_rateadj_cancer total_mortality_rateadj_diabendo total_mortality_rateadj_cvd total_mortality_rateadj_resp total_mortality_rateadj_digest total_mortality_rateadj_gu total_mortality_rateadj_psyc total_mortality_rateadj_external total_mortality_rateadj_other se_adj_matperi se_adj_infect se_adj_cancer se_adj_diabendo se_adj_cvd se_adj_resp se_adj_digest se_adj_gu se_adj_psyc se_adj_external se_adj_other, from(stdizemortsexbycause)

drop stdizemortsexbycause

label variable total_mortality_rate_adj_all "Total Mortality Rate (Adjusted to WHO Std) - Overall"
label variable total_mortality_rateadj_matperi "Total Mortality Rate (Adjusted to WHO Std) - Maternal & Perinatal Conditions"
label variable total_mortality_rateadj_infect "Total Mortality Rate (Adjusted to WHO Std) - Infectious Diseases"
label variable total_mortality_rateadj_cancer "Total Mortality Rate (Adjusted to WHO Std) - Cancer"
label variable total_mortality_rateadj_diabendo "Total Mortality Rate (Adjusted to WHO Std) - Diabetes & Endocrine Disorders"
label variable total_mortality_rateadj_cvd "Total Mortality Rate (Adjusted to WHO Std) - Circulatory Diseases"
label variable total_mortality_rateadj_resp "Total Mortality Rate (Adjusted to WHO Std) - Respiratory Diseases"
label variable total_mortality_rateadj_digest "Total Mortality Rate (Adjusted to WHO Std) - Digestive Diseases"
label variable total_mortality_rateadj_gu "Total Mortality Rate (Adjusted to WHO Std) - Genitourinary Diseases"
label variable total_mortality_rateadj_psyc "Total Mortality Rate (Adjusted to WHO Std) - Neuropsychiatric Conditions"
label variable total_mortality_rateadj_external "Total Mortality Rate (Adjusted to WHO Std) - External Causes"
label variable total_mortality_rateadj_other "Total Mortality Rate (Adjusted to WHO Std) - Other"

label variable total_ypll_rate_adj_all "Total YPLL Rate (Adjusted to WHO Std) - Overall"
label variable total_ypll_rateadj_matperi "Total YPLL Rate (Adjusted to WHO Std) - Maternal & Perinatal Conditions"
label variable total_ypll_rateadj_infect "Total YPLL Rate (Adjusted to WHO Std) - Infectious Diseases"
label variable total_ypll_rateadj_cancer "Total YPLL Rate (Adjusted to WHO Std) - Cancer"
label variable total_ypll_rateadj_diabendo "Total YPLL Rate (Adjusted to WHO Std) - Diabetes & Endocrine Disorders"
label variable total_ypll_rateadj_cvd "Total YPLL Rate (Adjusted to WHO Std) - Circulatory Diseases"
label variable total_ypll_rateadj_resp "Total YPLL Rate (Adjusted to WHO Std) - Respiratory Diseases"
label variable total_ypll_rateadj_digest "Total YPLL Rate (Adjusted to WHO Std) - Digestive Diseases"
label variable total_ypll_rateadj_gu "Total YPLL Rate (Adjusted to WHO Std) - Genitourinary Diseases"
label variable total_ypll_rateadj_psyc "Total YPLL Rate (Adjusted to WHO Std) - Neuropsychiatric Conditions"
label variable total_ypll_rateadj_external "Total YPLL Rate (Adjusted to WHO Std) - External Causes"
label variable total_ypll_rateadj_other "Total YPLL Rate (Adjusted to WHO Std) - Other"

label variable se_adj_all "Standard Error for Adjusted Mortality Rate - Overall"
label variable se_adj_matperi "Standard Error for Adjusted Mortality Rate - Maternal & Perinatal Conditions"
label variable se_adj_infect "Standard Error for Adjusted Mortality Rate - Infectious Diseases"
label variable se_adj_cancer "Standard Error for Adjusted Mortality Rate - Cancer"
label variable se_adj_diabendo "Standard Error for Adjusted Mortality Rate - Diabetes & Endocrine Disorders"
label variable se_adj_cvd "Standard Error for Adjusted Mortality Rate - Circulatory Diseases"
label variable se_adj_resp "Standard Error for Adjusted Mortality Rate - Respiratory Diseases"
label variable se_adj_digest "Standard Error for Adjusted Mortality Rate - Digestive Diseases"
label variable se_adj_gu "Standard Error for Adjusted Mortality Rate - Genitourinary Diseases"
label variable se_adj_psyc "Standard Error for Adjusted Mortality Rate - Neuropsychiatric Conditions"
label variable se_adj_external "Standard Error for Adjusted Mortality Rate - External Causes"
label variable se_adj_other "Standard Error for Adjusted Mortality Rate - Other"

label variable se_ypll_adj_all "Standard Error for Adjusted YPLL Rate - Overall"
label variable se_ypll_adj_matperi "Standard Error for Adjusted YPLL Rate - Maternal & Perinatal Conditions"
label variable se_ypll_adj_infect "Standard Error for Adjusted YPLL Rate - Infectious Diseases"
label variable se_ypll_adj_cancer "Standard Error for Adjusted YPLL Rate - Cancer"
label variable se_ypll_adj_diabendo "Standard Error for Adjusted YPLL Rate - Diabetes & Endocrine Disorders"
label variable se_ypll_adj_cvd "Standard Error for Adjusted YPLL Rate - Circulatory Diseases"
label variable se_ypll_adj_resp "Standard Error for Adjusted YPLL Rate - Respiratory Diseases"
label variable se_ypll_adj_digest "Standard Error for Adjusted YPLL Rate - Digestive Diseases"
label variable se_ypll_adj_gu "Standard Error for Adjusted YPLL Rate - Genitourinary Diseases"
label variable se_ypll_adj_psyc "Standard Error for Adjusted YPLL Rate - Neuropsychiatric Conditions"
label variable se_ypll_adj_external "Standard Error for Adjusted YPLL Rate - External Causes"
label variable se_ypll_adj_other "Standard Error for Adjusted YPLL Rate - Other"

order year sex total_mortality_rate_adj_all total_mortality_rateadj_matperi total_mortality_rateadj_infect total_mortality_rateadj_cancer  total_mortality_rateadj_diabendo total_mortality_rateadj_cvd total_mortality_rateadj_resp total_mortality_rateadj_digest total_mortality_rateadj_gu total_mortality_rateadj_psyc total_mortality_rateadj_external total_mortality_rateadj_other total_ypll_rate_adj_all total_ypll_rateadj_matperi total_ypll_rateadj_infect total_ypll_rateadj_cancer total_ypll_rateadj_diabendo total_ypll_rateadj_cvd total_ypll_rateadj_resp total_ypll_rateadj_digest total_ypll_rateadj_gu total_ypll_rateadj_psyc total_ypll_rateadj_external total_ypll_rateadj_other se_adj_all se_adj_matperi se_adj_infect se_adj_cancer se_adj_diabendo se_adj_cvd se_adj_resp se_adj_digest se_adj_gu se_adj_psyc se_adj_external se_adj_other se_ypll_adj_all se_ypll_adj_matperi se_ypll_adj_infect se_ypll_adj_cancer se_ypll_adj_diabendo se_ypll_adj_cvd se_ypll_adj_resp se_ypll_adj_digest se_ypll_adj_gu se_ypll_adj_psyc se_ypll_adj_external se_ypll_adj_other

save obj1_sex70_stdize, replace



/************************
 SECTION 3 - 
 BY AGE GROUP
************************/

/****SECTION 3.1 - BY AGE GROUP - OVERALL****/

frame create ageoverall
frame change ageoverall

use obj1_deaths_2010_2019

collapse ///
    (count) total_mortality_count = _temp ///
    (sum) total_ypll_count = ypll_70_burden ///
	(sum) total_ypll_squared = ypll_70_sq_burden, ///
    by(year agegroup)

drop if inlist(agegroup,7,8,.) /*no filling in needed -- all combos of age group and year represented*/

format total_mortality_count total_ypll_count %9.0gc

/*Add in population data*/

frame create populationdataage
frame change populationdataage

use "${clean_data_dir}/population_age.dta" 

frame change ageoverall

frlink 1:1 year age, frame(populationdataage)

frget population, from(populationdataage)


/*Calculate rates*/

generate total_mortality_rate = (total_mortality_count/population)*100000
generate total_ypll_rate = (total_ypll_count/population)*100000
generate se_ypll_rate = (sqrt(total_ypll_squared)/population)*100000

drop total_ypll_squared

format total_mortality_rate total_ypll_rate se_ypll_rate %9.0gc

order year agegroup population total_mortality_count total_mortality_rate total_ypll_count total_ypll_rate se_ypll_rate

rename (total_mortality_count total_mortality_rate total_ypll_count total_ypll_rate se_ypll_rate) (total_mortality_count_all total_mortality_rate_all total_ypll_count_all total_ypll_rate_all se_ypll_rate_all)

drop populationdataage


/****SECTION 3.2 - BY AGE - STRATIFIED BY MAJOR CAUSE****/

frame create agebymjrcause
frame change agebymjrcause

use obj1_deaths_2010_2019

collapse ///
    (count) total_mortality_count = _temp ///
    (sum) total_ypll_count = ypll_70_burden ///
    (sum) total_ypll_squared = ypll_70_sq_burden, ///
	by(year agegroup underlying_cause_cat)

drop if underlying_cause_cat==. | inlist(agegroup,7,8,.)

/*generate new rows for causes of death with 0 counts for a given age group - as these are missing from the dataset after the above collapse (due to 0 count)*/
fillin agegroup year underlying_cause_cat

replace total_mortality_count = 0 if missing(total_mortality_count)

replace total_ypll_count = 0 if missing(total_ypll_count)

replace total_ypll_squared = 0 if missing(total_ypll_squared)


format total_mortality_count total_ypll_count %9.0gc

drop _fillin

reshape wide total_mortality_count total_ypll_count total_ypll_squared, i(year agegroup) j(underlying_cause_cat) 

/*Add in population data and calculate rates*/

*overall population

frlink 1:1 year agegroup, frame(populationdataage)
frget population, from(populationdataage)


*calculate rates of mortality and premature mortality, rename and label variables

local codcat "matperi infect cancer diabendo cvd resp digest gu psyc external other"
local j = 1
forvalues j = 1/11 {
	local sffx : word `j' of `codcat' 
	rename (total_mortality_count`j' total_ypll_count`j' total_ypll_squared`j') (total_mortality_count_`sffx' total_ypll_count_`sffx' total_ypll_squared_`sffx')
    local ++j
}

local codcat2 "matperi infect cancer diabendo cvd resp digest gu psyc external other"
local q = 1
foreach x of local codcat2 {
    generate total_mortality_rate_`x' = (total_mortality_count_`x'/population)*100000
	generate total_ypll_rate_`x' = (total_ypll_count_`x'/population)*100000
    generate se_ypll_rate_`x' = (sqrt(total_ypll_squared_`x')/population)*100000
	format total_mortality_rate_`x' total_ypll_rate_`x' se_ypll_rate_`x' %9.0gc
    local ++q
}

drop total_ypll_squared_*

label variable total_mortality_count_matperi "Total Mortality Count - Maternal & Perinatal Conditions"
label variable total_mortality_rate_matperi "Total Mortality Rate - Maternal & Perinatal Conditions"
label variable total_ypll_count_matperi "Total YPLL Count - Maternal & Perinatal Conditions"
label variable total_ypll_rate_matperi "Total YPLL Rate - Maternal & Perinatal Conditions"
label variable se_ypll_rate_matperi "Standard Error for YPLL Rate - Maternal & Perinatal Conditions"

label variable total_mortality_count_infect "Total Mortality Count - Infectious Diseases"
label variable total_mortality_rate_infect "Total Mortality Rate - Infectious Diseases"
label variable total_ypll_count_infect "Total YPLL Count - Infectious Diseases"
label variable total_ypll_rate_infect "Total YPLL Rate - Infectious Diseases"
label variable se_ypll_rate_infect "Standard Error for YPLL Rate - Infectious Diseases"

label variable total_mortality_count_cancer "Total Mortality Count - Cancer"
label variable total_mortality_rate_cancer "Total Mortality Rate - Cancer"
label variable total_ypll_count_cancer "Total YPLL Count - Cancer"
label variable total_ypll_rate_cancer "Total YPLL Rate - Cancer"
label variable se_ypll_rate_cancer "Standard Error for YPLL Rate - Cancer"

label variable total_mortality_count_diabendo "Total Mortality Count - Diabetes & Endocrine Disorders"
label variable total_mortality_rate_diabendo "Total Mortality Rate - Diabetes & Endocrine Disorders"
label variable total_ypll_count_diabendo "Total YPLL Count - Diabetes & Endocrine Disorders"
label variable total_ypll_rate_diabendo "Total YPLL Rate - Diabetes & Endocrine Disorders"
label variable se_ypll_rate_diabendo "Standard Error for YPLL Rate - Diabetes & Endocrine Disorders"

label variable total_mortality_count_cvd "Total Mortality Count - Circulatory Diseases"
label variable total_mortality_rate_cvd "Total Mortality Rate - Circulatory Diseases"
label variable total_ypll_count_cvd "Total YPLL Count - Circulatory Diseases"
label variable total_ypll_rate_cvd "Total YPLL Rate - Circulatory Diseases"
label variable se_ypll_rate_cvd "Standard Error for YPLL Rate - Circulatory Diseases"

label variable total_mortality_count_resp "Total Mortality Count - Respiratory Diseases"
label variable total_mortality_rate_resp "Total Mortality Rate - Respiratory Diseases"
label variable total_ypll_count_resp "Total YPLL Count - Respiratory Diseases"
label variable total_ypll_rate_resp "Total YPLL Rate - Respiratory Diseases"
label variable se_ypll_rate_resp "Standard Error for YPLL Rate - Respiratory Diseases"

label variable total_mortality_count_digest "Total Mortality Count - Digestive Diseases"
label variable total_mortality_rate_digest "Total Mortality Rate - Digestive Diseases"
label variable total_ypll_count_digest "Total YPLL Count - Digestive Diseases"
label variable total_ypll_rate_digest "Total YPLL Rate - Digestive Diseases"
label variable se_ypll_rate_digest "Standard Error for YPLL Rate - Digestive Diseases"

label variable total_mortality_count_gu "Total Mortality Count - Genitourinary Diseases"
label variable total_mortality_rate_gu "Total Mortality Rate - Genitourinary Diseases"
label variable total_ypll_count_gu "Total YPLL Count - Genitourinary Diseases"
label variable total_ypll_rate_gu "Total YPLL Rate - Genitourinary Diseases"
label variable se_ypll_rate_gu "Standard Error for YPLL Rate - Genitourinary Diseases"

label variable total_mortality_count_psyc "Total Mortality Count - Neuropsychiatric Conditions"
label variable total_mortality_rate_psyc "Total Mortality Rate - Neuropsychiatric Conditions"
label variable total_ypll_count_psyc "Total YPLL Count - Neuropsychiatric Conditions"
label variable total_ypll_rate_psyc "Total YPLL Rate - Neuropsychiatric Conditions"
label variable se_ypll_rate_psyc "Standard Error for YPLL Rate - Neuropsychiatric Conditions"

label variable total_mortality_count_external "Total Mortality Count - External Causes"
label variable total_mortality_rate_external "Total Mortality Rate - External Causes"
label variable total_ypll_count_external "Total YPLL Count - External Causes"
label variable total_ypll_rate_external "Total YPLL Rate - External Causes"
label variable se_ypll_rate_external "Standard Error for YPLL Rate - External Causes"

label variable total_mortality_count_other "Total Mortality Count - Other"
label variable total_mortality_rate_other "Total Mortality Rate - Other"
label variable total_ypll_count_other "Total YPLL Count - Other"
label variable total_ypll_rate_other "Total YPLL Rate - Other"
label variable se_ypll_rate_other "Standard Error for YPLL Rate - Other"

frlink 1:1 year agegroup, frame(ageoverall)
frget total_mortality_count_all total_mortality_rate_all total_ypll_count_all total_ypll_rate_all se_ypll_rate_all, from(ageoverall)

label variable total_mortality_count_all "Total Mortality Count - Overall"
label variable total_mortality_rate_all "Total Mortality Rate - Overall"
label variable total_ypll_count_all "Total YPLL Count - Overall"
label variable total_ypll_rate_all "Total YPLL Rate - Overall"
label variable se_ypll_rate_all "Standard Error for YPLL Rate - Overall"

order year agegroup population total_mortality_count_all total_mortality_rate_all total_mortality_count_matperi total_mortality_rate_matperi total_mortality_count_infect total_mortality_rate_infect total_mortality_count_cancer total_mortality_rate_cancer total_mortality_count_diabendo total_mortality_rate_diabendo total_mortality_count_cvd total_mortality_rate_cvd total_mortality_count_resp total_mortality_rate_resp total_mortality_count_digest total_mortality_rate_digest total_mortality_count_gu total_mortality_rate_gu total_mortality_count_psyc total_mortality_rate_psyc total_mortality_count_external total_mortality_rate_external total_mortality_count_other total_mortality_rate_other total_ypll_count_all total_ypll_rate_all total_ypll_count_matperi total_ypll_rate_matperi total_ypll_count_infect total_ypll_rate_infect total_ypll_count_cancer total_ypll_rate_cancer total_ypll_count_diabendo total_ypll_rate_diabendo total_ypll_count_cvd total_ypll_rate_cvd total_ypll_count_resp total_ypll_rate_resp total_ypll_count_digest total_ypll_rate_digest total_ypll_count_gu total_ypll_rate_gu total_ypll_count_psyc total_ypll_rate_psyc total_ypll_count_external total_ypll_rate_external total_ypll_count_other total_ypll_rate_other se_ypll_rate_all se_ypll_rate_matperi se_ypll_rate_infect se_ypll_rate_cancer se_ypll_rate_diabendo se_ypll_rate_cvd se_ypll_rate_resp se_ypll_rate_digest se_ypll_rate_gu se_ypll_rate_psyc se_ypll_rate_external se_ypll_rate_other  

sort year agegroup

drop populationdataage ageoverall 

save obj1_age70, replace


/************************
 SECTION 4 - 
  BY AGE GROUP & SEX
************************/

/****SECTION 4.1 - BY AGE GROUP & SEX - OVERALL****/

frame create agesexoverall
frame change agesexoverall

use obj1_deaths_2010_2019

collapse ///
    (count) total_mortality_count = _temp ///
    (sum) total_ypll_count = ypll_70_burden ///
	(sum) total_ypll_squared = ypll_70_sq_burden, ///
    by(year agegroup sex)

drop if inlist(agegroup,7,8,.) | sex==9 /*no filling in needed -- all combos of age group and year represented*/

format total_mortality_count total_ypll_count %9.0gc

/*Add in population data*/

frame create populationdataagesex
frame change populationdataagesex

use "${clean_data_dir}/population_agesex.dta" 

frame change agesexoverall

frlink 1:1 year age sex, frame(populationdataagesex)

frget population, from(populationdataagesex)


/*Calculate rates*/

generate total_mortality_rate = (total_mortality_count/population)*100000
generate total_ypll_rate = (total_ypll_count/population)*100000
generate se_ypll_rate = (sqrt(total_ypll_squared)/population)*100000

drop total_ypll_squared

format total_mortality_rate total_ypll_rate se_ypll_rate %9.0gc

order year agegroup sex population total_mortality_count total_mortality_rate total_ypll_count total_ypll_rate se_ypll_rate

rename (total_mortality_count total_mortality_rate total_ypll_count total_ypll_rate se_ypll_rate) (total_mortality_count_all total_mortality_rate_all total_ypll_count_all total_ypll_rate_all se_ypll_rate_all)

drop populationdataagesex


/*SECTION 4.2 - BY AGE GROUP & SEX - STRATIFIED BY MAJOR CAUSE */

frame create agesexbymjrcause
frame change agesexbymjrcause

use obj1_deaths_2010_2019

collapse ///
    (count) total_mortality_count = _temp ///
    (sum) total_ypll_count = ypll_70_burden ///
	(sum) total_ypll_squared = ypll_70_sq_burden, ///
    by(year agegroup sex underlying_cause_cat)

drop if underlying_cause_cat==. | inlist(agegroup,7,8,.) | sex==9


/*generate new rows for causes of death with 0 counts for a given age group - as these are missing from the dataset after the above collapse (due to 0 count)*/
fillin agegroup sex year underlying_cause_cat

replace total_mortality_count = 0 if missing(total_mortality_count)

replace total_ypll_count = 0 if missing(total_ypll_count)

replace total_ypll_squared = 0 if missing(total_ypll_squared)


format total_mortality_count total_ypll_count %9.0gc

drop _fillin

reshape wide total_mortality_count total_ypll_count total_ypll_squared, i(year agegroup sex) j(underlying_cause_cat) 

/*Add in population data and calculate rates*/

frlink 1:1 year agegroup sex, frame(populationdataagesex)
frget population, from(populationdataagesex)


*calculate rates of mortality and premature mortality, rename and label variables

local codcat "matperi infect cancer diabendo cvd resp digest gu psyc external other"
local j = 1
forvalues j = 1/11 {
	local sffx : word `j' of `codcat' 
	rename (total_mortality_count`j' total_ypll_count`j' total_ypll_squared`j') (total_mortality_count_`sffx' total_ypll_count_`sffx' total_ypll_squared_`sffx')
    local ++j
}

local codcat2 "matperi infect cancer diabendo cvd resp digest gu psyc external other"
local q = 1
foreach x of local codcat2 {
    generate total_mortality_rate_`x' = (total_mortality_count_`x'/population)*100000
	generate total_ypll_rate_`x' = (total_ypll_count_`x'/population)*100000
	generate se_ypll_rate_`x' = (sqrt(total_ypll_squared_`x')/population)*100000
	format total_mortality_rate_`x' total_ypll_rate_`x' se_ypll_rate_`x' %9.0gc
    local ++q
}

drop total_ypll_squared_*

label variable total_mortality_count_matperi "Total Mortality Count - Maternal & Perinatal Conditions"
label variable total_mortality_rate_matperi "Total Mortality Rate - Maternal & Perinatal Conditions"
label variable total_ypll_count_matperi "Total YPLL Count - Maternal & Perinatal Conditions"
label variable total_ypll_rate_matperi "Total YPLL Rate - Maternal & Perinatal Conditions"
label variable se_ypll_rate_matperi "Standard Error for YPLL Rate - Maternal & Perinatal Conditions"

label variable total_mortality_count_infect "Total Mortality Count - Infectious Diseases"
label variable total_mortality_rate_infect "Total Mortality Rate - Infectious Diseases"
label variable total_ypll_count_infect "Total YPLL Count - Infectious Diseases"
label variable total_ypll_rate_infect "Total YPLL Rate - Infectious Diseases"
label variable se_ypll_rate_infect "Standard Error for YPLL Rate - Infectious Diseases"

label variable total_mortality_count_cancer "Total Mortality Count - Cancer"
label variable total_mortality_rate_cancer "Total Mortality Rate - Cancer"
label variable total_ypll_count_cancer "Total YPLL Count - Cancer"
label variable total_ypll_rate_cancer "Total YPLL Rate - Cancer"
label variable se_ypll_rate_cancer "Standard Error for YPLL Rate - Cancer"

label variable total_mortality_count_diabendo "Total Mortality Count - Diabetes & Endocrine Disorders"
label variable total_mortality_rate_diabendo "Total Mortality Rate - Diabetes & Endocrine Disorders"
label variable total_ypll_count_diabendo "Total YPLL Count - Diabetes & Endocrine Disorders"
label variable total_ypll_rate_diabendo "Total YPLL Rate - Diabetes & Endocrine Disorders"
label variable se_ypll_rate_diabendo "Standard Error for YPLL Rate - Diabetes & Endocrine Disorders"

label variable total_mortality_count_cvd "Total Mortality Count - Circulatory Diseases"
label variable total_mortality_rate_cvd "Total Mortality Rate - Circulatory Diseases"
label variable total_ypll_count_cvd "Total YPLL Count - Circulatory Diseases"
label variable total_ypll_rate_cvd "Total YPLL Rate - Circulatory Diseases"
label variable se_ypll_rate_cvd "Standard Error for YPLL Rate - Circulatory Diseases"

label variable total_mortality_count_resp "Total Mortality Count - Respiratory Diseases"
label variable total_mortality_rate_resp "Total Mortality Rate - Respiratory Diseases"
label variable total_ypll_count_resp "Total YPLL Count - Respiratory Diseases"
label variable total_ypll_rate_resp "Total YPLL Rate - Respiratory Diseases"
label variable se_ypll_rate_resp "Standard Error for YPLL Rate - Respiratory Diseases"

label variable total_mortality_count_digest "Total Mortality Count - Digestive Diseases"
label variable total_mortality_rate_digest "Total Mortality Rate - Digestive Diseases"
label variable total_ypll_count_digest "Total YPLL Count - Digestive Diseases"
label variable total_ypll_rate_digest "Total YPLL Rate - Digestive Diseases"
label variable se_ypll_rate_digest "Standard Error for YPLL Rate - Digestive Diseases"

label variable total_mortality_count_gu "Total Mortality Count - Genitourinary Diseases"
label variable total_mortality_rate_gu "Total Mortality Rate - Genitourinary Diseases"
label variable total_ypll_count_gu "Total YPLL Count - Genitourinary Diseases"
label variable total_ypll_rate_gu "Total YPLL Rate - Genitourinary Diseases"
label variable se_ypll_rate_gu "Standard Error for YPLL Rate - Genitourinary Diseases"

label variable total_mortality_count_psyc "Total Mortality Count - Neuropsychiatric Conditions"
label variable total_mortality_rate_psyc "Total Mortality Rate - Neuropsychiatric Conditions"
label variable total_ypll_count_psyc "Total YPLL Count - Neuropsychiatric Conditions"
label variable total_ypll_rate_psyc "Total YPLL Rate - Neuropsychiatric Conditions"
label variable se_ypll_rate_psyc "Standard Error for YPLL Rate - Neuropsychiatric Conditions"

label variable total_mortality_count_external "Total Mortality Count - External Causes"
label variable total_mortality_rate_external "Total Mortality Rate - External Causes"
label variable total_ypll_count_external "Total YPLL Count - External Causes"
label variable total_ypll_rate_external "Total YPLL Rate - External Causes"
label variable se_ypll_rate_external "Standard Error for YPLL Rate - External Causes"

label variable total_mortality_count_other "Total Mortality Count - Other"
label variable total_mortality_rate_other "Total Mortality Rate - Other"
label variable total_ypll_count_other "Total YPLL Count - Other"
label variable total_ypll_rate_other "Total YPLL Rate - Other"
label variable se_ypll_rate_other "Standard Error for YPLL Rate - Other"

frlink 1:1 year agegroup sex, frame(agesexoverall)
frget total_mortality_count_all total_mortality_rate_all total_ypll_count_all total_ypll_rate_all se_ypll_rate_all, from(agesexoverall)

label variable total_mortality_count_all "Total Mortality Count - Overall"
label variable total_mortality_rate_all "Total Mortality Rate - Overall"
label variable total_ypll_count_all "Total YPLL Count - Overall"
label variable total_ypll_rate_all "Total YPLL Rate - Overall"
label variable se_ypll_rate_all "Standard Error for YPLL Rate - Overall"

order year agegroup sex population total_mortality_count_all total_mortality_rate_all total_mortality_count_matperi total_mortality_rate_matperi total_mortality_count_infect total_mortality_rate_infect total_mortality_count_cancer total_mortality_rate_cancer total_mortality_count_diabendo total_mortality_rate_diabendo total_mortality_count_cvd total_mortality_rate_cvd total_mortality_count_resp total_mortality_rate_resp total_mortality_count_digest total_mortality_rate_digest total_mortality_count_gu total_mortality_rate_gu total_mortality_count_psyc total_mortality_rate_psyc total_mortality_count_external total_mortality_rate_external total_mortality_count_other total_mortality_rate_other total_ypll_count_all total_ypll_rate_all total_ypll_count_matperi total_ypll_rate_matperi total_ypll_count_infect total_ypll_rate_infect total_ypll_count_cancer total_ypll_rate_cancer total_ypll_count_diabendo total_ypll_rate_diabendo total_ypll_count_cvd total_ypll_rate_cvd total_ypll_count_resp total_ypll_rate_resp total_ypll_count_digest total_ypll_rate_digest total_ypll_count_gu total_ypll_rate_gu total_ypll_count_psyc total_ypll_rate_psyc total_ypll_count_external total_ypll_rate_external total_ypll_count_other total_ypll_rate_other se_ypll_rate_all se_ypll_rate_matperi se_ypll_rate_infect se_ypll_rate_cancer se_ypll_rate_diabendo se_ypll_rate_cvd se_ypll_rate_resp se_ypll_rate_digest se_ypll_rate_gu se_ypll_rate_psyc se_ypll_rate_external se_ypll_rate_other  

sort year agegroup sex

drop populationdataagesex agesexoverall

save obj1_agesex70, replace

/************************
 SECTION 5 - 
 ASSESSING DIFFERENCES 
 IN MEDIAN YPLL AND
 PCT OF PREMATURE DEATHS
 BY SEX FOR 2010 AND 2019
************************/
use obj1_deaths_2010_2019
describe

/* 1  - ASSESS DIFFERENCES IN PERCENT OF PREMATURE DEATHS BY SEX AND CAUSE  */

// Step 1: Decode cause labels if not already done
decode underlying_cause_cat, gen(cause_str)

// Step 2: Define years of interest
local years 2010 2019

// Step 3: Create postfile with percent premature added
tempname memhold
tempfile results_all
postfile `memhold' str40 cause year chi2_stat p_value n_obs pct_premature_f pct_premature_m premature_female premature_male using `results_all'

// Step 4: Loop over years
foreach y of local years {

    levelsof cause_str if year == `y', local(causes)

    // --- All causes ---
    quietly count if year == `y' & sex == 1
    local Nf = r(N)
    quietly count if year == `y' & sex == 1 & premature70 == 1
    local Pf = r(N)
    local pctf = 100 * `Pf' / `Nf'

    quietly count if year == `y' & sex == 2
    local Nm = r(N)
    quietly count if year == `y' & sex == 2 & premature70 == 1
    local Pm = r(N)
    local pctm = 100 * `Pm' / `Nm'

    quietly count if year == `y' & inlist(sex, 1, 2)
    local n = r(N)

    quietly tabulate sex premature70 if year == `y' & inlist(sex, 1, 2), chi2
    if r(p) < . {
        post `memhold' ("All causes") (`y') (r(chi2)) (r(p)) (`n') (`pctf') (`pctm') (`Pf') (`Pm')
    }

    // --- By cause ---
    foreach c of local causes {

        quietly count if year == `y' & cause_str == "`c'" & sex == 1
        local Nf = r(N)
        quietly count if year == `y' & cause_str == "`c'" & sex == 1 & premature70 == 1
        local Pf = r(N)
        local pctf = cond(`Nf' > 0, 100 * `Pf' / `Nf', .)

        quietly count if year == `y' & cause_str == "`c'" & sex == 2
        local Nm = r(N)
        quietly count if year == `y' & cause_str == "`c'" & sex == 2 & premature70 == 1
        local Pm = r(N)
        local pctm = cond(`Nm' > 0, 100 * `Pm' / `Nm', .)

        quietly count if year == `y' & cause_str == "`c'" & inlist(sex, 1, 2)
        local n = r(N)

        quietly tabulate sex premature70 if year == `y' & cause_str == "`c'" & inlist(sex, 1, 2), chi2
        if r(p) < . {
            post `memhold' ("`c'") (`y') (r(chi2)) (r(p)) (`n') (`pctf') (`pctm') (`Pf') (`Pm')
        }
    }
}

// Step 5: Finalize and save
postclose `memhold'
use `results_all', clear

// Step 6: Bonferroni correction
gen n_tests = _N
gen p_bonf = p_value * n_tests
replace p_bonf = 1 if p_bonf > 1

// Step 7: Save
save obj1_chi2results_70, replace



/* 2  - ASSESS DIFFERENCES IN MEAN YPLL BY SEX AND CAUSE */

use obj1_deaths_2010_2019, clear
describe

/* histogram ypll_70 if sex == 1, normal title("YPLL for Females")
histogram ypll_70 if sex == 2, normal title("YPLL for Males")*/

summarize ypll_70 if sex == 1, detail
summarize ypll_70 if sex == 2, detail  

/*******************************************************
NOTE AUGUST 1, 2025 - THE YPLL VALUES ARE NOT 
NORMALLY DISTRIBUTED WITHIN EACH SEX. 
THEREFORE T-TEST CANNOT BE RUN 
AND IT MAY NOT BE APPROPRIATE TO REPORT THE MEAN.
WILL HAVE TO ADD MEDIAN YPLL TO CALCULATIONS ABOVE 
AND RUN THE MANN-WHITNEY U / WILCOXON RANK-SUM TEST
******************************************************/

// ASSESS DIFFERENCES IN MEDIAN YPLL

// --- Step 1. Decode cause labels if needed ---
decode underlying_cause_cat, gen(cause_str)

// --- Step 2. Define years you want to analyze ---
local years 2010 2019

// --- Step 3. Prepare postfile ---
tempname memhold
tempfile results_all2

postfile `memhold' str40 cause year median_female median_male N_female N_male p_value using `results_all2'

// --- Step 4. Loop over years ---
foreach z of local years {

    // --- Get list of causes for this year ---
    levelsof cause_str if year == `z', local(causes2)

    // --- Loop over causes ---
    foreach c of local causes2 {
        
        // Get counts
        quietly count if year == `z' & cause_str == "`c'" & sex == 1 & !missing(ypll_70)
        local N_female = r(N)
        
        quietly count if year == `z' & cause_str == "`c'" & sex == 2 & !missing(ypll_70)
        local N_male = r(N)

        // Skip if either group has <2
        if `N_male' < 2 | `N_female' < 2 {
            continue
        }

        // Compute medians
        centile ypll_70 if year == `z' & cause_str == "`c'" & sex == 1, centile(50)
        local med_female = r(c_1)

        centile ypll_70 if year == `z' & cause_str == "`c'" & sex == 2, centile(50)
        local med_male = r(c_1)

        // Wilcoxon rank-sum
        quietly ranksum ypll_70 if year == `z' & cause_str == "`c'" & inlist(sex,1,2), by(sex)
        local p = r(p)

        // Post results
        post `memhold' ("`c'") (`z') (`med_female') (`med_male') (`N_female') (`N_male') (`p')
    }

    // --- Add overall (all causes) for this year ---
    quietly count if year == `z' & sex == 1 & !missing(ypll_70)
    local N_female = r(N)

    quietly count if year == `z' & sex == 2 & !missing(ypll_70)
    local N_male = r(N)

    if `N_male' >= 2 & `N_female' >= 2 {
        centile ypll_70 if year == `z' & sex == 1, centile(50)
        local med_female = r(c_1)

        centile ypll_70 if year == `z' & sex == 2, centile(50)
        local med_male = r(c_1)

        quietly ranksum ypll_70 if year == `z' & inlist(sex,1,2), by(sex)
        local p = r(p)

        post `memhold' ("All causes") (`z') (`med_female') (`med_male') (`N_female') (`N_male') (`p')
    }
}

postclose `memhold'

// --- Step 5. Load results and compute Bonferroni-adjusted p-values ---
use `results_all2', clear
gen n_tests = _N
gen p_bonf = p_value * n_tests
replace p_bonf = 1 if p_bonf > 1

// --- Step 6. Save results ---
save obj1_wmwresults_70, replace


/*export datasets to excel*/

local exportlist obj1_overall70 obj1_overall70_stdize obj1_sex70 obj1_sex70_stdize obj1_age70 obj1_agesex70 obj1_chi2results_70 obj1_wmwresults_70
foreach ds of local exportlist {
    use `ds', clear
    export excel using "${analysis_output_dir}/Output - Age 70.xlsx", ///
    sheet("`ds'", replace) firstrow(variables)

/*prepare the above datasets for joinpoint regression in the NCI joinpoint software

first, reshape each dataset from wide to long. include covariates (age, sex, cause of death), independent variable (year), and the dependent variable (YPLL rate). 
next, calcualte the standard error for rates for which this has not yet been done above (i.e., the crude and age-specific rates)

then, append:
1. overall and by sex (crude) datasets
2. by age dataset
3. overall and by sex (age standardized) datasets 

*/

**overall - crude**
frame create dataforjoinpointall
frame change dataforjoinpointall

use obj1_overall70

keep year population70 total_ypll_count_all total_ypll_rate_all se_ypll_rate_all total_ypll_count_matperi total_ypll_rate_matperi se_ypll_rate_matperi total_ypll_count_infect total_ypll_rate_infect se_ypll_rate_infect total_ypll_count_cancer total_ypll_rate_cancer se_ypll_rate_cancer total_ypll_count_diabendo total_ypll_rate_diabendo se_ypll_rate_diabendo total_ypll_count_cvd total_ypll_rate_cvd se_ypll_rate_cvd total_ypll_count_resp total_ypll_rate_resp se_ypll_rate_resp total_ypll_count_digest total_ypll_rate_digest se_ypll_rate_digest total_ypll_count_gu total_ypll_rate_gu se_ypll_rate_gu total_ypll_count_psyc total_ypll_rate_psyc se_ypll_rate_psyc total_ypll_count_external total_ypll_rate_external se_ypll_rate_external total_ypll_count_other total_ypll_rate_other se_ypll_rate_other            

reshape long total_ypll_count_@ total_ypll_rate_@ se_ypll_rate_@, i(year) j(cause) string

gen sex="ALL"

replace cause=strupper(cause)

rename (total_ypll_count_ total_ypll_rate_ se_ypll_rate_) (ypll_count ypll_rate se_ypll_rate)

save obj1_overall70_reshape, replace


**sex - crude**
frame create dataforjoinpointsex
frame change dataforjoinpointsex

use obj1_sex70

keep year population70 sex total_ypll_count_all total_ypll_rate_all se_ypll_rate_all total_ypll_count_matperi total_ypll_rate_matperi se_ypll_rate_matperi total_ypll_count_infect total_ypll_rate_infect se_ypll_rate_infect total_ypll_count_cancer total_ypll_rate_cancer se_ypll_rate_cancer total_ypll_count_diabendo total_ypll_rate_diabendo se_ypll_rate_diabendo total_ypll_count_cvd total_ypll_rate_cvd se_ypll_rate_cvd total_ypll_count_resp total_ypll_rate_resp se_ypll_rate_resp total_ypll_count_digest total_ypll_rate_digest se_ypll_rate_digest total_ypll_count_gu total_ypll_rate_gu se_ypll_rate_gu total_ypll_count_psyc total_ypll_rate_psyc se_ypll_rate_psyc total_ypll_count_external total_ypll_rate_external se_ypll_rate_external total_ypll_count_other total_ypll_rate_other se_ypll_rate_other

reshape long total_ypll_count_@ total_ypll_rate_@ se_ypll_rate_@, i(year sex) j(cause) string

decode sex, generate(sex2)

drop sex

replace sex2=strupper(sex2)

rename sex2 sex

replace cause=strupper(cause)

rename (total_ypll_count_ total_ypll_rate_ se_ypll_rate_) (ypll_count ypll_rate se_ypll_rate)

order cause sex year ypll_count ypll_rate se_ypll_rate
sort cause sex year

save obj1_sex70_reshape, replace

append using obj1_overall70_reshape, keep(cause sex year ypll_count ypll_rate se_ypll_rate population70)

sort cause sex year
order cause sex year population70 ypll_count ypll_rate se_ypll_rate 

save obj1_overallsex70_crude_jnpt, replace

export delimited using obj1_overallsex70_crude_jnpt, replace


**age**
frame create dataforjoinpointage
frame change dataforjoinpointage

use obj1_age70

keep year population agegroup total_ypll_count_all total_ypll_rate_all se_ypll_rate_all total_ypll_count_matperi total_ypll_rate_matperi se_ypll_rate_matperi total_ypll_count_infect total_ypll_rate_infect se_ypll_rate_infect total_ypll_count_cancer total_ypll_rate_cancer se_ypll_rate_cancer total_ypll_count_diabendo total_ypll_rate_diabendo se_ypll_rate_diabendo total_ypll_count_cvd total_ypll_rate_cvd se_ypll_rate_cvd total_ypll_count_resp total_ypll_rate_resp se_ypll_rate_resp total_ypll_count_digest total_ypll_rate_digest se_ypll_rate_digest total_ypll_count_gu total_ypll_rate_gu se_ypll_rate_gu total_ypll_count_psyc total_ypll_rate_psyc se_ypll_rate_psyc total_ypll_count_external total_ypll_rate_external se_ypll_rate_external total_ypll_count_other total_ypll_rate_other se_ypll_rate_other

reshape long total_ypll_count_@ total_ypll_rate_@ se_ypll_rate_@, i(year agegroup) j(cause) string

decode agegroup, generate(agegroup2)

drop agegroup

replace agegroup2=strupper(agegroup2)

rename agegroup2 agegroup

replace cause=strupper(cause)

rename (total_ypll_count_ total_ypll_rate_ se_ypll_rate_) (ypll_count ypll_rate se_ypll_rate)

order cause agegroup year population ypll_count ypll_rate se_ypll_rate
sort cause agegroup year

save obj1_age70_crude_jnpt, replace

export delimited using obj1_age70_crude_jnpt, replace


**overall - standardized**
frame create dataforjoinpointallstd
frame change dataforjoinpointallstd

use obj1_overall70_stdize

keep year total_ypll_rate_adj_all total_ypll_rateadj_matperi total_ypll_rateadj_infect total_ypll_rateadj_cancer total_ypll_rateadj_diabendo total_ypll_rateadj_cvd total_ypll_rateadj_resp total_ypll_rateadj_digest total_ypll_rateadj_gu total_ypll_rateadj_psyc total_ypll_rateadj_external total_ypll_rateadj_other se_ypll_adj_all se_ypll_adj_matperi se_ypll_adj_infect se_ypll_adj_cancer se_ypll_adj_diabendo se_ypll_adj_cvd se_ypll_adj_resp se_ypll_adj_digest se_ypll_adj_gu se_ypll_adj_psyc se_ypll_adj_external se_ypll_adj_other

rename total_ypll_rate_adj_all total_ypll_rateadj_all

reshape long total_ypll_rateadj_@ se_ypll_adj_@, i(year) j(cause) string

gen sex="ALL"

replace cause=strupper(cause)
generate cause2=subinstr(cause, "ADJ", "", .)

drop cause

rename cause2 cause

rename (total_ypll_rate se_ypll_adj_) (ypll_rate_adj se_ypll_adj)

order cause sex year ypll_rate_adj se_ypll_adj
sort cause sex year

save obj1_overall70_stdize_reshape, replace


**sex - standardized**
frame create dataforjoinpointsexstd
frame change dataforjoinpointsexstd

use obj1_sex70_stdize

keep year sex total_ypll_rate_adj_all total_ypll_rateadj_matperi total_ypll_rateadj_infect total_ypll_rateadj_cancer total_ypll_rateadj_diabendo total_ypll_rateadj_cvd total_ypll_rateadj_resp total_ypll_rateadj_digest total_ypll_rateadj_gu total_ypll_rateadj_psyc total_ypll_rateadj_external total_ypll_rateadj_other se_ypll_adj_all se_ypll_adj_matperi se_ypll_adj_infect se_ypll_adj_cancer se_ypll_adj_diabendo se_ypll_adj_cvd se_ypll_adj_resp se_ypll_adj_digest se_ypll_adj_gu se_ypll_adj_psyc se_ypll_adj_external se_ypll_adj_other

rename total_ypll_rate_adj_all total_ypll_rateadj_all

reshape long total_ypll_rateadj_@ se_ypll_adj_@, i(year sex) j(cause) string

decode sex, generate(sex2)
drop sex
replace sex2=strupper(sex2)

rename sex2 sex

replace cause=strupper(cause)
generate cause2=subinstr(cause, "ADJ", "", .)
drop cause
rename cause2 cause

rename (total_ypll_rate se_ypll_adj_) (ypll_rate_adj se_ypll_adj)

order cause sex year ypll_rate_adj se_ypll_adj
sort cause sex year

save obj1_sex70_stdize_reshape, replace

append using obj1_overall70_stdize_reshape, keep(cause sex year ypll_rate_adj se_ypll_adj)

sort cause sex year

order cause sex year ypll_rate_adj se_ypll_adj

save obj1_overallsex70_stdize_jnpt, replace

export delimited using obj1_overallsex70_stdize_jnpt, replace





/*******************************************************************************
RGD vs JCF Homicide Comparison
Purpose: Cross-validate RGD-recorded homicide deaths against JCF-reported 
         homicide victim counts, 2013-2019, ages 0-64 years
*******************************************************************************/

use obj1_deaths_2010_2019, clear

/*restrict to homicide ICD-10 codes (X85-Y09) and study years 2013-2019*/
keep if inrange(year, 2013, 2019)

tab1 age, m
keep if inrange(age, 0, 64) /*ages 0-64 only to align with JCF data)*/
tab1 age, m

* 1. Standardize by stripping out periods, spaces, and forcing uppercase
gen clean_icd = upper( ///
    subinstr( ///
    subinstr(trim(underlying_cause_clean), ".", "", .), ///
    " ", "", .))

* 2. Generate strict and expanded versions of homicides
gen byte homicide_strict = regexm(clean_icd, ///
    "^(X8[5-9]|X9[0-9]|Y0[0-9]|Y871)")

gen byte homicide_expanded = homicide_strict | ///
    regexm(clean_icd, "^(Y22|Y23|Y24|Y28|Y29)")

assert homicide_expanded >= homicide_strict

tab1 homicide_strict homicide_expanded, m
tab1 underlying_cause_clean if homicide_strict==1
tab1 underlying_cause_clean if homicide_expanded==1
tab1 underlying_cause_clean, m

* 3. Generate age groups specific to this analysis
recode age ///
    (0/14  = 1) ///
    (15/29 = 2) ///
    (30/39 = 3) ///
    (40/49 = 4) ///
    (50/59 = 5) ///
    (60/64 = 6), gen(jcf_agegroup)
	
label define jcf_agegroup ///
    1 "0 to 14" ///
    2 "15 to 29" ///
    3 "30 to 39" ///
    4 "40 to 49" ///
    5 "50 to 59" ///
    6 "60 to 64"

label values jcf_agegroup jcf_agegroup
assert !missing(jcf_agegroup)

* 4. Annual counts 
preserve
collapse (sum) ///
    rgd_strict=homicide_strict ///
    rgd_expanded=homicide_expanded, by(year)

isid year
assert rgd_expanded >= rgd_strict
save rgd_homicide_counts_by_year, replace
restore

* 5. Age-specific counts 
collapse (sum) ///
    rgd_strict=homicide_strict ///
    rgd_expanded=homicide_expanded, by(jcf_agegroup)

isid jcf_agegroup
assert rgd_expanded >= rgd_strict
save rgd_homicide_counts_by_age, replace

/*******************************************************************************
JCF Homicide Trends
Purpose: Examine changes in trends in homicides reported by JCF from 2013-2019,
		ages 0-64 years
*******************************************************************************/
capture frame drop jcftrends
frame create jcftrends
frame change jcftrends

* Import JCF data
import excel "${input_dir}/jcf_murders_2013_2019.xlsx", ///
    sheet("Sheet1") firstrow case(lower)

describe
assert inrange(year, 2013, 2019)
assert inrange(agegroup, 1, 13)
assert !missing(year, agegroup, murders)
assert murders >= 0

isid year agegroup
assert _N == 91
bysort year: assert _N == 13

*Estimate YPLL based on the midpoint of each five-year age group
gen age_midpoint = (agegroup*5) - 2.5
gen ypll_per_death = 70 - age_midpoint
gen jcf_ypll = murders*ypll_per_death

rename agegroup agegroup_5yr
save jcf_murders_2013_2019, replace


*Prepare population data 
use "${clean_data_dir}/population_age_5yrbands.dta", clear

tab1 year agegroup_5yr, m 
tab1 agegroup_5yr, m nolabel

keep if inrange(year, 2013, 2019) & inrange(agegroup_5yr, 1, 13)

assert !missing(year, agegroup_5yr, population)
assert population > 0

isid year agegroup_5yr
assert _N == 91
bysort year: assert _N == 13

*Merge and calculate rate
merge 1:1 year agegroup_5yr using jcf_murders_2013_2019 

tab1 _merge
assert _merge == 3
drop _merge


/*merge in WHO standard population weights (for subset of pop <65 years)*/
capture frame drop weights65
frame create weights65
frame change weights65

use "${clean_data_dir}/whostandard_65.dta" 

isid agegroup_5yr
assert _N == 13
assert inrange(agegroup_5yr, 1, 13)
assert !missing(weight_whostd_65)
assert weight_whostd_65 > 0

frame change jcftrends

frlink m:1 agegroup_5yr, frame(weights65)

frget weight_whostd_65, from(weights65)

drop weights65

tab1 agegroup_5yr, nolabel 

/*calculate standardized rates*/
generate crude = jcf_ypll / population
generate product = crude*weight_whostd_65
gen double variance_squared = (weight_whostd_65^2 * murders * (ypll_per_death^2)) / (population^2)


/*sum across age groups*/
collapse ///
    (sum) _jcf_ypll_rate_adj = product ///
    (sum) _se_ypll_adj= variance_squared, ///
    by(year)

generate jcf_ypll_rate_adj = _jcf_ypll_rate_adj * 100000
gen se_ypll_adj = (sqrt(_se_ypll_adj)) * 100000

format jcf_ypll_rate_adj se_ypll_adj %12.4f

drop _jcf_ypll_rate_adj _se_ypll_adj


*Save and export
save jcf_murder_trend_2013_2019, replace
export delimited using "jcf_murder_joinpoint_2013_2019.csv", replace

/*******************************************************************************
SECTION 6
SUPPLEMENTARY POISSON REGRESSION ANALYSES

Purpose:
Estimate age- and year-adjusted differences in:
1. YPLL rates by sex
2. Premature mortality rates by sex

Unit of analysis:
Year x 5-year age group x sex
*******************************************************************************/

/**************************************************
STEP 1. CREATE POPULATION DENOMINATOR MASTER
**************************************************/

frame create poisson_models
frame change poisson_models

use "${clean_data_dir}/population_agesex_5yrbands.dta" 

tab1 agegroup_5yr sex, m

keep if inrange(agegroup_5yr, 1, 14)

isid year agegroup_5yr sex

* Check denominator
assert population > 0
assert !missing(population)

tempfile population_poisson
save `population_poisson'

/**************************************************
STEP 2. CREATE NUMERATORS
**************************************************/

use obj1_deaths_2010_2019, clear

* Restrict to analytic strata
keep if inrange(agegroup_5yr, 1, 14)
keep if inlist(sex, 1, 2)

collapse ///
    (sum) total_ypll = ypll_70_burden ///
    (sum) premature_deaths = premature70, ///
    by(year agegroup_5yr sex)
	
isid year agegroup_5yr sex

tempfile deaths_poisson
save `deaths_poisson'

/**************************************************
STEP 3. MERGE NUMERATORS INTO POPULATION MASTER
**************************************************/

use `population_poisson', clear

merge 1:1 year agegroup_5yr sex using `deaths_poisson'

tab _merge

* There should be no death strata without a population denominator
assert _merge != 2
assert _N == 280

* Population strata with no events are real zero-event strata
replace total_ypll = 0 if _merge == 1
replace premature_deaths = 0 if _merge == 1

drop _merge

assert !missing(total_ypll)
assert !missing(premature_deaths)

isid year agegroup_5yr sex

save obj1_poisson_models_70, replace


/* MODEL 1: YPLL */
use obj1_poisson_models_70, clear

poisson total_ypll ib1.sex i.agegroup_5yr i.year, ///
    exposure(population) ///
    vce(robust) ///
    irr

estimates store ypll_model

test 2.sex
estat gof, pearson


/* MODEL 2: PREMATURE DEATHS */
poisson premature_deaths ib1.sex i.agegroup_5yr i.year, ///
    exposure(population) ///
    vce(robust) ///
    irr

estimates store premature_model

test 2.sex
estat gof, pearson



/*******************************************************************************
Descriptive Table 1
Characteristics of registered deaths and premature mortality,
Jamaica, 2010–2019
*******************************************************************************/
use obj1_deaths_2010_2019, clear
preserve


/**************************SEX VARIABLE FOR COLUMNS**************************/
generate byte table1_sex = .
replace table1_sex = 1 if sex == 1
replace table1_sex = 2 if sex == 2
replace table1_sex = 3 if !inlist(sex, 1, 2)

label define table1_sexlbl ///
    1 "Female" ///
    2 "Male" ///
    3 "Missing/unknown", replace

label values table1_sex table1_sexlbl
label variable table1_sex "Sex"

tab sex table1_sex, missing


/* Premature death count and percentage */
clonevar table1_premature = premature70

label define table1_premlbl 0 "No" 1 "Yes", replace
label values table1_premature table1_premlbl
label variable table1_premature "Deaths before age 70, n (%)"


/* Separate copies allow total and median YPLL to appear on different rows */
clonevar table1_ypll_total = ypll_70
clonevar table1_ypll_median = ypll_70

label variable table1_ypll_total "Total YPLL"
label variable table1_ypll_median ///
    "Median YPLL per premature death (IQR)"

/* Premature death age groups*/

generate byte table1_agegroup = .

replace table1_agegroup = 1 if premature70 == 1 & inrange(age,  0, 14)
replace table1_agegroup = 2 if premature70 == 1 & inrange(age, 15, 29)
replace table1_agegroup = 3 if premature70 == 1 & inrange(age, 30, 39)
replace table1_agegroup = 4 if premature70 == 1 & inrange(age, 40, 49)
replace table1_agegroup = 5 if premature70 == 1 & inrange(age, 50, 59)
replace table1_agegroup = 6 if premature70 == 1 & inrange(age, 60, 69)

label define table1_agelbl ///
    1 "0–14" ///
    2 "15–29" ///
    3 "30–39" ///
    4 "40–49" ///
    5 "50–59" ///
    6 "60–69", replace

label values table1_agegroup table1_agelbl
label variable table1_agegroup ///
    "Age at death among premature deaths, n (%)"

assert !missing(table1_agegroup) if premature70 == 1
assert missing(table1_agegroup) if premature70 == 0

generate double table1_age_premature = age if premature70 == 1

label variable table1_age_premature ///
    "Age at premature death, median (IQR), years"

assert missing(table1_age_premature) if premature70 == 0
assert table1_age_premature == age if premature70 == 1

/* Causes of premature deaths*/

generate byte table1_cause = underlying_cause_cat if premature70 == 1

label define table1_causelbl ///
     1 "Maternal and perinatal conditions" ///
     2 "Infectious and parasitic diseases" ///
     3 "Cancer" ///
     4 "Diabetes and endocrine disorders" ///
     5 "Circulatory diseases" ///
     6 "Chronic respiratory diseases" ///
     7 "Oral and digestive diseases" ///
     8 "Genitourinary diseases" ///
     9 "Neuropsychiatric conditions" ///
    10 "External causes" ///
    11 "Other causes", replace

label values table1_cause table1_causelbl
label variable table1_cause ///
    "Underlying cause among premature deaths, n (%)"

assert !missing(table1_cause) if premature70 == 1
assert missing(table1_cause) if premature70 == 0

tab table1_cause if premature70 == 1, missing


collect clear

dtable ///
    age ///
    1.table1_premature ///
    table1_ypll_total ///
    table1_ypll_median ///
	table1_age_premature ///
    i.table1_agegroup ///
    i.table1_cause, ///
    by(table1_sex, notests totals) ///
    sample("Registered deaths, n", ///
        statistics(frequency) place(items)) ///
    continuous(age, ///
        statistics(q2 iqi)) ///
    continuous(table1_ypll_total, ///
        statistics(total)) ///
    continuous(table1_ypll_median, ///
        statistics(q2 iqi)) ///
	continuous(table1_age_premature, ///
    statistics(q2 iqi)) ///
    factor(table1_premature table1_agegroup table1_cause, ///
        statistics(fvfrequency fvpercent) test(none)) ///
    define(iqi = q1 q3, delimiter("–")) ///
    nformat(%12.0fc frequency fvfrequency total) ///
    nformat(%9.1f mean sd fvpercent q1 q2 q3 iqi) ///
    sformat("(%s)" sd) ///
    sformat("(%s%%)" fvpercent) ///
    sformat("[%s]" iqi) ///
    column(by(hide) total("Overall")) ///
    title("Table 1. Characteristics of registered deaths and premature mortality in Jamaica, 2010–2019.") ///
	note ("Values represent N (%) unless otherwise specified.") ///
    note("Female and male columns exclude records with missing sex; these records remain in the overall column.") ///
    note("YPLL = years of potential life lost before age 70; IQR = interquartile range.") ///
	export("${manuscript_output_dir}/Table1.xlsx", replace)
	
	
capture log close _all
	




