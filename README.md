## Analysis scripts used for manuscript: 

###Title: Resource acquisition and allocation traits in symbiotic rhizobia with implications for life-history outside of legume hosts

#### Authors:
Katherine E. Muller (collected data)
R. Ford Denison


--------------------------
Script:
## PHBvariationInNaturalPopulations_edited.R

Requirements:
packages: plyr, ggplot2, lme4, sjstats
data files:
###chamaecristaCombinedNodData_edited.RDATA
### senescentPlusJulyNodPHB_Becker.RDATA
### SS3regatedEditedPHBdata.RData
### phbplusplantdata_final.RDATA
### editedPHBData_WREsplitroot1.RDATA
### noduleMassWithH2assays.RDATA
### isolateCodesForFreezerwithfilenames_WRE.csv
### flowDataA3B5nod_bdpy.csv
### PHBstandards_A3B5.csv
### 20150206biomassA3B5.csv

Sections:
###Results (main text)
#### Phenotypic distribution of PHB accumulation in natural rhizobia populations
#### Heritable phenotypic variation in PHB accumulation
### code for Figure 1: PHB accumulation varies widely among co-occurring rhizobia
### Appendix A: Additional details on materials and methods
#### Statistical tests for heritability estimates
### Appendix C: Additional results on genetic an environmental factors contributing to variation in rhizobial PHB accumulation. 
### Appendix D: Among- and within-strain PHB variation measured in singly inoculated soybean plants. 
------------------------
Script:
## simulatingPHBuse_edited.R
Requirements: 
scripts: 
###functionsForEstimatingSurvivalOnPHB.R

data: 
###senescentNoduleDataEditedA3B1.RData

Sections:
## Main text
### Results
#### Model predictions

## Appendix B
 
---------------------------------------
Script:
## soilTempGraphs.R
Requirements:
data: averagemonthlysoiltemp_waseca.csv

Sections:
### Appendix B
#### Figure B1: Seasonal oscillation in soil temperature

---------------------------------------
Script:
## starvingNoduleExtracts_final.R

Requirements:
packages: ggplot2
scripts:
###functionsForEstimatingSurvivalOnPHB.R

data files:
### starvingNoduleExtractPHBuse.RDATA
### viabilityCountsEdited.RData

Sections:

## Results (main text)
### Population trends in starving rhizobia
### PHB use in starving rhizobia
### Figure 3: Initial PHB only partially explains population trends in starvation culture.
### Figure 4: Segregation of high- and low-PHB subgroups in starvation culture suggest asymmetric partitioning of PHB to old-pole cells.
### Figure 5: PHB use during starvation

### Appendix E

-------------------------------------------
Script:
## starvingNoduleExtract_contaminants.R

Requirements:
data:
### starvingNoduleExtractPHBuse.RData
### contaminationNST1.csv

Sections:
### Appendix E
