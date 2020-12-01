The current folder contains files for implementing the four-level GEE approach introduced in “Power Considerations for GEE Analyses of Four-Level Clustered Designs” by X. Wang et al.

Binary outcomes are generated from a marginal mean model with specified correlation structures to mimic a four-level cluster randomized trial. The intervention is assigned at the cluster level.

List of Files:
1) binGEN.R = R file for simulating binary outcome data in a balanced four-level CRT
2) binGEN_var.R = R file for simulating binary outcome data in an unbalanced four-level CRT
3) binMAEE.R = GEE and MAEE program for clustered binary outcomes
4) LogitBinBCV.R = GEE program for clustered binary outcomes
5) example_analysis.R = R file that illustrates the use of binMAEE/LogitBinBCV main program based on a simulated dataset using binGEN/binGEN_var

NOTES: 1) The example program demonstrates the computation for a trial with 14 clusters, 2 divisions per cluster, 3 subjects per division, and 5 evaluations per subject (or on average). 2) You will need to change path/directory names to import the functions before running the example program.
