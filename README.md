The current folder includes R code for reproducing all of the tables and figures in the article “Power considerations for generalized estimating equations analyses of four-level cluster randomized trials” by Wang et al. and in the Supporting Information for this article. It also includes a file of the Supporting Information for this article.

For questions or comments about the code please contact Xueqi Wang at xueqi.wang@duke.edu.

I. Supporting Information of the Main Article:

1. Supporting information for four-level CRTs, by Wang et al.pdf = file that includes the Web Appendices, Tables, and Figures referenced in the main article.

II. List of Supporting Files: These supporting files are sourced in the main files that reproduce the numbers in the submitted manuscript and supporting information.

2. binGEN.R = function to simulate binary outcome data in a balanced four-level CRT;
3. binGEN_var.R = function to simulate binary outcome data in an unbalanced four-level CRT;
4. binMAEE.R = function to calculate the bias-corrected sandwich variances for logistic-binomial GEE/MAEE analyses;
5. LogitBinBCV.R = function to calculate the bias-corrected sandwich variances for logistic-binomial GEE analyses with the independence working correlation model.

III. List of Main Files: These main files are used to reproduce the results in the submitted manuscript and supporting information.

6. binScenarios_Analytical.R = reproduce predicted power results in all of the tables and figures for power;
7. binScenarios_MAEE_00.R = reproduce simulation results of GEE/MAEE analyses using the extended nested exchangeable working correlation structure, under balanced four-level CRTs (used in Tables 2-3, Figures 2-3, Web Table 1, and Web Figures 1-8);
8. binScenarios_MAEE_25.R = reproduce simulation results of GEE/MAEE analyses using the extended nested exchangeable working correlation structure, under unbalanced four-level CRTs with CV = 0.25 (used in Figures 2-3, Web Tables 4-5, and Web Figures 1-2);
9. binScenarios_MAEE_50.R = reproduce simulation results of GEE/MAEE analyses using the extended nested exchangeable working correlation structure, under unbalanced four-level CRTs with CV = 0.50 (used in Figures 2-3, Web Tables 6-7, and Web Figures 3-4);
10. binScenarios_MAEE_75.R = reproduce simulation results of GEE/MAEE analyses using the extended nested exchangeable working correlation structure, under unbalanced four-level CRTs with CV = 0.75 (used in Figures 2-3, Web Tables 8-9, and Web Figures 5-6);
11. binScenarios_MAEE_100.R = reproduce simulation results of GEE/MAEE analyses using the extended nested exchangeable working correlation structure, under unbalanced four-level CRTs with CV = 1.00 (used in Figures 2-3, Web Tables 10-11, and Web Figures 7-8);
12. binScenarios_ind.R = reproduce simulation results of GEE analyses using an independence working correlation matrix, under both balanced and unbalanced four-level CRTs (used in Figures 2-3, Web Tables 2-3, Web Tables 12-19, and Web Figures 1-8);
13. binRES.R = reproduce numbers and the format of Tables 2-3, Figures 2-3, Web Tables 1-19, and Web Figures 1-8, after running the above simulation programs;
14. Application.R = reproduce results of applications in Figure 4-5, Web Figures 9-12, and Web Tables 20-21.

IV. Folder

15. binResults = folder to save power/size data (from the raw data through binScenarios_Analytical.R and binRES.R) with an subfolder “binRData,” which contains 30 empty sub-subfolders ("1row" to "30row") to save raw data from simulation programs 7-12.

V. Software

Analyses were conducted with R, version 4.0.2 (https://www.r-project.org/). The calculations used R packages MASS (version 7.3-51.6), gee (version 4.13-20), openxlsx (version 4.1.5), ggplot2 (version 3.3.2), directlabels (version 2020.6.17), and cowplot (version 1.0.0).

VI. R commands for the installation of R packages

install.packages(c("MASS", "gee", "openxlsx", "ggplot2", "directlabels", "cowplot"))

NOTES: Make sure the current working directory is the current folder before running each program; the folders mentioned in #15 should be created before running the programs.
