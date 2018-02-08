# euc_rates

### This repository contains all the necessary information for reproducing the results of "Faster rates of evolution in taller eucalypts with larger genomes". 

This paper looked at four correlates of molecular rate variation in eucalypts: specific leaf area (SLA), height, seed mass and genome size.

The raw data files are in the "data" folder. This includes:

1. Sequence alignments 
- a. Non-coding nuclear sequences 1391 base pairs long, for 711 species
- b. Non-coding chloroplast sequences 1767 base pairs long, for 711 species
- c. whole chloroplast genome sequences 139598 base pairs long, for 41 species

2. Trait data for SLA, height, seed mass and genome size (with accepted species names checked agains the Australian Plant Census). Included as separate files are two novel datasets:
- a. SLA measurements for 511 eucalypt species
- b. genome size measurements for 788 eucalypt species

3. phylogenetic trees created with IQ-TREE from the sequence alignments

This analysis was performed in both a ML and Bayesian framework and all of the necessary scripts are available.

### ML analysis
- Input: phylogenetic trees and trait data
- Output: summary statistics of correlation between substitution rates and species traits

Method:
1. Tree is pruned to remove branches with fewer than 5 substitutions as these are known to make large differences to substitution rate estimates
2. Time-calibrated phylogeny created (both with and without fossils)
3. Substitution rates are extracted from the terminal branches of the time-calibrated phylogeny
4. Substitution rate and trait data is combined into a single data.frame
5. PGLS analysis is used to look for correlations between substitution rates and traits

To use our data files, please change the path in our script to point to where you have saved our files.

To reproduce our results:
1. The prune.tree() function will be used throughout the analysis to prune branch lengths with less than 5 substitutions. Use the 'prune_short_branches_fn' file to create these functions in your global environment to use whenever they are required in the script.
2. Begin with the 'cross_validation_lambda' file to determine appropriate values for lambda. (Lambda is the rate smoothing parameter. It is used when creating a time calibrated phylogeny)
3. Use the 'ML_analysis' file to determine correlations between rates and traits without fossil data for the trees 'cp_tree.phy' and 'nuc_tree.phy'. This file randomly prunes ML trees 200 times to avoid biasing correlations between rates and traits. Repeat this analysis for both the chloroplast tree and the nuclear tree.
4. Use the 'ML_analysis_fossils' file to determine correlations between rates and traits with fossil data for the trees 'cp_tree.phy' and 'nuc_tree.phy'. This file randomly prunes ML trees 200 times to avoid biasing correlations between rates and traits. Repeat this analysis for both the chloroplast tree and the nuclear tree.
5. Use the 'ML_wholecp' file to determine correlations between rates and traits without fossil data for the 'whole_cptree.tre' file. This script can also be used to calculate correlations between averaged traits by pointing at the 'LHS_averaged' file instead of the 'merged_LHS' file.
6. Use the 'ML_wholecp_fossils' file to determine correlations between rates and traits with fossil data for the 'whole_cptree.tre' file. This script can also be used to calculate correlations between averaged traits by pointing at the 'LHS_averaged' file instead of the 'merged_LHS' file. 
7. For each of these analyses, save your 'rates' files to a separate folder as you will use them later.
8. Use the 'ML_rate_variation' file to visually inspect the distribution of rate estimates. 

Because the rate estimates we obtained did not agree with those found elsewhere in the literature, we repeated the analysis in a Bayesian framework.

### Bayesian analysis (coevol)
- Input: phylogenetic tree topology, sequence alignments, trait data
- Output: table of correlations, covariances and posterior probabilities between substitution rates and traits

Method:
1. To reduce computational time, files were pruned to contain less than 100 species
2. Save each alignment file as a CSV file
3. The prune.tree() function will be used throughout the analysis to prune branch lengths with less than 5 substitutions. Use the 'prune_short_branches_fn' file to create these functions in your global environment to use whenever they are required in the script.
4. Use the 'coevol_setup' file to:
- a. Remove branches with <5 substitutions
- b. Remove species for which complete trait data is not available
- c. Continue to prune out species until 100 remain
- Then, use the species in the tree file to prune unwanted species from the alignment and trait files. Repeat this for all three datasets.
- The 'LHS_average' file is used only for the whole chloroplast data
5. Drop mallees from the analysis using the 'coevol_setup_mallees' file. This is similar to the previous step except that mallee species are also excluded.
6. Alignment files then need to be manually changed back into phylip format. Remove all " symbols.
7. For analyses involving fossils, the 'fossils' file needs to be altered so that so that nodes are specified by by giving the names of two species that have this node as their last common ancestor. This can be done by visually inspecting the tree file.
8. Use these files to run the program coevol. The commandline prompt is in the 'commandline.sh' file. Separate files need to be created for each analysis. That is, for each dataset: with and without fossils and with and without mallees. 
