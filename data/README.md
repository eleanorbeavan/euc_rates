# This file explains the contents of each file in the /data folder

#### 1. LHS_averaged.csv

- This file contains life-history trait data for each species in the whole chloroplast dataset averaged over that lineage in the four-gene dataset

  - **Column 1:** accepted_name. The name of the species given by the Australian Plant Census https://biodiversity.org.au/nsl/services/apc
  - **Column 2:** coevol_name. The name of each species with the genus name truncated for the Coevol program
  - **Column 3:** mean_gs. The lineage-averaged genome size measurement (2C value) for each species 
  - **Column 4:** mean_height. The lineage-averaged height measurement (m) for each species
  - **Column 5:** mean_sla. The lineage-averaged SLA measurement (cm2/g) for each species
  - **Column 6:** mean_seed. The lineage-averaged seed mass measurement (g) for each species

#### 2. alignment_names.csv

- This file provides a reference for species names in the alignment, trait and coevol files

  - **Column 1:** alignment_name. The name of each species in the two-gene alignment files
  - **Column 2:** accepted_name. The name of each species given by the Australian Plant Census https://biodiversity.org.au/nsl/services/apc
  - **Column 3:** coevol. The name of each species with the genus name truncated for the Coevol program

#### 3. cp_alignment.nex

- A nexus file with the two-gene chloroplast alignment

#### 4. cp_tree.phy

- A phylip tree file with the two-gene chloroplast phylogeny obtained from IQ-TREE using the cp_alignment.nex file

#### 5. four_gene_tree.tre

- A tree file with the four-gene nuclear and chloroplast phylogeny obtained from IQ-TREE using the concatenated nuc_alignment.nex and cp_alignment.nex files

#### 6. genome_size.csv

- This file contains the raw genome size measurement data 

  - **Column 1:** accepted_name. The name of each species given by the Australian Plant Census https://biodiversity.org.au/nsl/services/apc
  - **Column 2:** species_name. The name of the species as given by Currency Creek Arboretum
  - **Column 3:** Voucher specimen. The number assigned to the sample
  - **Column 4:** Generalised location collected. The location where the seed from the mother tree was collected in the native range of that species.
  - **Column 5:** Row. The row number of the tree at Currency Creek Arboretum
  - **Column 6:** Tree. The tree number of the tree at Currency Creek Arboretum
  - **Column 7:** Year planted. The year the tree was planted
  - **Column 8:** 1st budded. The year the tree first produced buds
  - **Column 9:** 1st flowered. The year the tree first produced flowers
  - **Column 10:** date collected. The date the samples were collected
  - **Column 11:** index. The genome size value recorded for the reference species Solanum pseudocapsicum used to calibrate the genome size measurements.
  - **Column 12:** genome_size. The genome size measurement for each species
  - **Column 13:** CV sample. Coefficient of variation values for the samples
  - **Column 14:** CV Solanum. Coefficient of variation values for the reference species Solanum pseudocapsicum
  - **Column 15:** fcm date. The data the flow cell cytometry was performed. 
  
#### 7. mallee.csv

- This file provides information on which species are mallees based on the Euclid database

  - **Column 1:** species_name. The name of each species given by Euclid
  - **Column 2:** accepted_name. The name of each species given by the Australian Plant Census https://biodiversity.org.au/nsl/services/apc
  - **Column 3:** coevol. The name of each species with the genus name truncated for the Coevol program
  - **Column 4:** mallee. The habit of each species. N = never mallee, S = sometimes mallee, Y = always mallee.


#### 8. merged_LHS.csv

- This file contains life-history trait data (genome size, height, SLA, and seed mass) for all available eucalypt species. Where multiple trait values were availabe for a species in the raw data sets the average was taken. If multiple values were attributed to more than one subspecies, we searched the literature for those subspecies and used the trait from the most common subspecies if one subspecies was far more common than the others, or took the average of the trait values of all subspecies if all available subspecies were roughly equally represented in the present-day population. 

  - **Column 1:** accepted_name. The name of the species given by the Australian Plant Census https://biodiversity.org.au/nsl/services/apc
  - **Column 2:** coevol_name. The name of each species with the genus name truncated for the Coevol program
  - **Column 3:** mean_gs. The final genome size measurement (2C value) for each species 
  - **Column 4:** mean_height. The final height measurement (m) for each species
  - **Column 5:** mean_sla. The final SLA measurement (cm2/g) for each species
  - **Column 6:** mean_sm. The final seed mass measurement (g) for each species
  - **Column 7:** genus. The genus of each species

#### 9. nuc_alignment.nex

- A nexus file with the two-gene nuclear alignment

#### 10. nuc_tree.phy

- A phylip tree file with the two-gene nuclear phylogeny obtained from IQ-TREE using the nuc_alignment.nex file

#### 11. sla_eucalypt.csv

- This file contains the raw SLA measurement data 

  - **Column 1:** Sample (Bag). The sample number
  - **Column 2:** accepted_name. The name of the species given by the Australian Plant Census https://biodiversity.org.au/nsl/services/apc
  - **Column 3:** species. The name of the species as given by Currency Creek Arboretum
  - **Column 4-6:** dry mass (g). The dry mass of each leaf replicate (1,2,3) after being oven dried at 45 degrees for 72 hours.
  - **Column 7-9:** Area mm2. The area of each leaf replicate (1,2,3) measured using ImageJ.
  - **Column 10-12:** SLA cm2/g. The SLA value for each leaf replicate (1,2,3). Calculated by converting the leaf area to cm2, by dividing by 100 and then dividing this area by leaf mass for each replicate.
  - **Column 13:** mean_SLA. The mean of the three SLA values for each sample.
  - **Column 14:** CV SLA
  - **Column 15:** Mean Mass. The average value for mass for the three replicates
  - **Column 16:** Mean Area. The average value for area for the three replicates
 
  
#### 12. whole_cptree.tre

- A tree file with the whole chloroplast phylogeny obtained from IQ-TREE using the wholecp_alignment.nex file

#### 13. wholecp_alignment.nex

- A nexus file with the whole chloroplast alignment

