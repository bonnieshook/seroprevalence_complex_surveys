/******************************************************************************
/* Program: example_code.sas
/* Creation Date: 01 FEB 2023
/* Manuscript: "Accounting for outcome misclassification in complex sample 
	designs: estimation of SARS-CoV-2 seroprevalence in central North Carolina"
/* Purpose: This SAS file provides code to calculate corrected seroprevalence
	estimates with corresponding 95% bootstrap confidence intervals				
/* Authors: Caitlin Cassidy, Bonnie Shook-Sa
/* Inputs: sample file with components listed in corresponding README file
*******************************************************************************/

*** Specify estimated sensitivity and specificity of diagnostic test;
%let sens=0.897; 
%let spec=0.993; 

*** Specify the sample sizes used to estimate sensitivity and specificity;
%let n_sens=145;
%let n_spec=274;

*** Specify number of bootstrap samples;
%let nboot=100; 

*** Replace with file path of example data set;
%let path = /example_data.csv;

*** Read in hypothetical cohort of 200 individuals with prevalence approx. 20% - 
see README document for more information;
*** Replace with own dataset of similar structure;
proc import datafile = "&path" 
	out=example
	dbms=csv
	replace;
run;

title "First 15 observations of dataset";
proc print data = example (obs=15);
run;

*** Calibrate weights, calculate unadjusted estimate (p_hat_s - Eqn 1);
proc sort data=example; 
	by STRATUM CLUSTER_NUM; 
run;
*** See SUDAAN Language Manual for more information about this procedure;
*** Population counts in POSTWGT statement come from American Community Survey 
	estimates (2020);
PROC WTADJUST DATA=example DESIGN=WR ADJUST=POST MAXITER=200 P_EPSILON=.001;
	WEIGHT BASEWT;
	NEST STRATUM CLUSTER_NUM;
	CLASS age_cat sex;
 	MODEL _one_=age_cat sex;
 	/*overall target population size then marginals for age group and sex in order 
 	of PROC FREQ output*/
 	POSTWGT 60970 25042 35928 32121 28849; 
 	VAR serostatus;
 	IDVAR participant_id;
 	OUTPUT MEAN SE_MEAN LOW_MEAN UP_MEAN / filename=uncorrected_est replace;
 	OUTPUT / predicted=all filename=calibrated_wts filetype=sas replace; 
run;

*Merge calibrated weights with example dataset;
proc sort data=example; 
	by participant_id; 
run;
proc sort data=calibrated_wts; 
	by participant_id; 
run;
data example2;
 	merge example(in=a)
       calibrated_wts(in=b keep=participant_id WTFINAL);
	by participant_id;
run;

*confirm that weights sum to marginal population totals used in 
POSTWGT statement in WTADJUST procedure;
title "";
proc means data=example2 min median max sum; 
	class age_cat; 
	var WTFINAL; 
run;
proc means data=example2 min median max sum; 
	class sex; 
	var WTFINAL; 
run;

*** Calculate Rogan-Gladen estimate (p_hat_c - Eqn 2);
data uncorrected_est2;
	set uncorrected_est(rename=(mean=mean_test) where=(_C1=0));
	uncor_prev=mean_test;
	uncor_prev_SE = se_mean;
	uncor_prev_LL=Low_mean;
	uncor_prev_UL=up_mean;
	cor_prev = (uncor_prev+&spec.-1)/(&sens.+&spec.-1); *calculate p' (R-G adjusted estimate);
	if cor_prev < 0 then cor_prev = 0; *truncate to 0;
	if cor_prev > 1 then cor_prev = 1; *truncate to 1;
	call symput('phat_c', cor_prev);
	keep uncor_prev uncor_prev_SE cor_prev uncor_prev_LL uncor_prev_UL ; *keep p and p' and CIs for uncorrected p;
run;

*** Calculate stratum sizes for bootstrap samples;
proc sort data=example2 out=clusters nodupkey; 
	by CLUSTER_NUM STRATUM; 
run;
proc freq data=clusters noprint; 
	tables STRATUM / list out=clusters2(drop=percent); 
run;

*** Select bootstrap samples for the computation of confidence intervals;
*We will sample m_h-1 PSUs with replacement from each stratum;
data strat_samsize_boot;
  	set clusters2;
 	_nsize_=(count-1)*&nboot.; 
run;

proc sort data=example2 out=unique_clusters(keep=CLUSTER_NUM STRATUM PSU_WT PER_WT) nodupkey; 
	by CLUSTER_NUM STRATUM; 
run;
proc sort data=unique_clusters; 
	by STRATUM;
run;

***Rao-Wu rescaling bootstrap;
*select bootstrap samples with replacement;
proc surveyselect data=unique_clusters method=urs n=strat_samsize_boot out = boot_smps noprint; 
	strata STRATUM;
run;

*** Divide into random bootstrap samples;
data boot_smps2;
 	set boot_smps;
 	do i=1 to NumberHits;
 		random=ranuni(12345);
 		output;
 	end;
run;

proc sort data=boot_smps2; 
	by STRATUM random; 
run;

data boot_smps3;
 	set boot_smps2;
 	by STRATUM;
 	if first.STRATUM then BOOT_NUM=0;
 	if BOOT_NUM GE &nboot. then BOOT_NUM=0;
 	BOOT_NUM+1;
run;

*Collapse file to the PSU level within bootstrap samples;
proc freq data=boot_smps3 noprint; 
	tables BOOT_NUM*STRATUM*CLUSTER_NUM*PSU_WT / list out=boot_smps4(drop=percent rename=(count=NUM_HITS)); 
run;

*define bootstrap weights;
proc sort data=boot_smps4; 
	by STRATUM CLUSTER_NUM; 
run;

data boot_smps5(rename=(CLUSTER_NUM=PSU));
 	merge boot_smps4(in=a) clusters2(in=b);
 	by STRATUM;
 	boot_weight=PSU_WT*(count/(count-1))*NUM_HITS;
run;

proc sort data=boot_smps5; 
	by BOOT_NUM; 
run;

*merge back with PERWT and HHWT;
proc sql;
 	create table boot_smps5b as
 	select A.*, B.participant_id, B.serostatus, B.HH_WT, B.PER_WT, B.CLUSTER_NUM, B.sex, B.age_cat
 	from boot_smps5 A left join example2 B
 	on A.PSU=B.CLUSTER_NUM 
 	order by CLUSTER_NUM;
quit;

proc sort data= boot_smps5b; 
	by BOOT_NUM STRATUM CLUSTER_NUM participant_id;
run;

data boot_smps6;
 	set boot_smps5b;
 	by BOOT_NUM;
 	if first.BOOT_NUM then CASECOUNT=0;
 	CASECOUNT+1;
 	FINWT_boot = boot_weight*HH_WT*PER_WT*15; *multiply by 15 to get in range of chatham pop size pre-calibration;
run;


*** Calibrate bootstrap weights for each replicate sample;
%macro calibrateboot;
%do iter=1 %to &nboot.;
data iter&iter.;
 	set boot_smps6;
 	if BOOT_NUM=&iter.;
run;

*Calibrate bootstrap weights for each individual (Eqn 3);
PROC WTADJUST DATA=iter&iter. DESIGN=WR ADJUST=POST MAXITER=500 P_EPSILON=.001;
 	WEIGHT FINWT_boot;
 	NEST _one_;
 	CLASS age_cat sex;
 	MODEL _one_=age_cat sex;
 	POSTWGT 60970 25042 35928 32121 28849; /*overall pop size then marginals for age group and sex*/
 	VAR serostatus;
 	IDVAR CASECOUNT;
 	OUTPUT MEAN / filename=uncorrected_est_iter&iter. replace; 
run;

*Calculate p_hat_s based on each replicate sample b (Eqn 4);
data uncorrected2_est_iter&iter.;
 	set uncorrected_est_iter&iter.;
 	if _C1=0;
 	that=mean;
 	iter=&iter.;
 	keep iter that;
run;
%end;
%mend calibrateboot;

%calibrateboot;

*** Combine and compute p_hat_c for each replicate b (Eqn 5);
*** Randomly draw sensitivity and specificity for each bootstrap sample;
data boot_ests;
 	set uncorrected2_est_iter1-uncorrected2_est_iter&nboot.;
 	sens_iter_boot=Rand("BINOMIAL", &sens., &n_sens.)/&n_sens.;
 	spec_iter_boot=Rand("BINOMIAL", &spec., &n_spec.)/&n_spec.;
 	phat=(that+spec_iter_boot-1)/(sens_iter_boot+spec_iter_boot-1);
 	if phat<0 then phat=0;
 	if phat>1 then phat=1;
run;

*** Calculate bootstrap CI as the 2.5 and 97.5th quantiles of the bootstrap estimates;
proc univariate data=boot_ests noprint;
  	var phat;
  	output out=boot_cis pctlpre=P_ pctlpts= 2.5 97.5;
run;

*store limits as macro variables;
data boot_cis2;
	set boot_cis;
	call symput('low', P_2_5);
	call symput('high', P_97_5);
run;

*print results to results window;
title "Corrected Prevalence with 95% Bootstrap CI";
data results;
	CorrectedPrevalence = &phat_c;
	LowerCLBoot = &low;
	UpperCLBoot = &high;
run;

proc print data = results;
run;
