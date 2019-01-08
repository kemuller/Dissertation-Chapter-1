## Analysis scripts used for manuscript:

### Title: Resource acquisition and allocation traits in symbiotic rhizobia with implications for life-history outside of legume hosts

#### Authors:
Katherine E. Muller (collected data)
R. Ford Denison

#### Date of data collection:

#### Geographic location of data collection (where was data collected?):

#### Information about funding sources that supported the collection of the data:


--------------------------
Script:
## PHBvariationInNaturalPopulations_edited.R

Requirements:
packages: plyr, ggplot2, lme4, sjstats
data files:
### chamaecristaCombinedNodData_edited.RDATA
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
### Results (main text)
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
### functionsForEstimatingSurvivalOnPHB.R

data:
### senescentNoduleDataEditedA3B1.RData

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

-----------------------------------------

### File specific information for datasets

#### PHBstandards_A3B5.csv
	file
	FLmean
	gatedCells
	date
	stain
	standardID
	PHBperCellpg
	standardID2
	PHBperCellpg2

#### averagemonthlysoiltemp_waseca.csv
	year
	month
	average4inchtempF
	average8inchtempF
	average4inchtempC
	average8inchtempC

#### contaminationNST1.csv
	sample
	date
	density
	strain
	type
	count1
	r1
	count2
	r2
	count3
	r3
	v1
	v2
	v3
	rep

#### flowDataA3B5nod_bdpy.csv
	file
	FL1mean
	gatedCells
	date
	plantID
	nodAgeDays
	nodLet
	dilution
	note
	delete

#### isolateCodesForFreezerwithfilenames_WRE.csv
	nodID
	file
	phbEst
	code
	repFromGrowthCurves
	plantID
	treatment
	rotation
	phbEstfromArithmeticMean

#### chamecristaCombinedNodData_edited.RDAtA:
##### nod.data.combined
	plantID
	nodID
	file
	FL1mean
	FL1cv
	FSCmean
	FSCcv
	SSCmean
	SSCcv
	countGated
	FL1mean_high
	FSCmean_high
	SSCmean_high
	FL1mean_low
	FSCmean_low
	SSCmean_low
	countGated_high
	countGated_low
	runDate
	dilutionFromNodule
	dilutionFromTube
	flowRate
	note
	phbEst
	wt.g
	cfuperml
	propSmooth
	growsWithBG
	growthStage
	shootDWg

#### editedPHBData_WREsplitroot1.RDATA
##### phbdata
	file
	FL1geomean
	FL1mode
	FL1cv
	FSCmean
	FSCmode
	FSCcv
	SSCmean
	SSCmode
	SSCcv
	countGated
	countUngated
	FSCgeomean
	FL1mean
	X
	runDate
	plate
	row
	column
	strain
	plantID
	skip
	dilution
	flowRate_ulpermin
	note
	phbEst


##### standards
	sampleID
	file
	FL1geomean
	FL1mode
	FL1cv
	FSCmean
	FSCmode
	FSCcv
	SSCmean
	SSCmode
	SSCcv
	countGated
	countUngated
	FSCgeomean
	FL1mean
	X
	runDate
	plate
	row
	column
	plantID
	skip
	dilution
	flowRate_ulpermin
	note
	phbPerCellPg


#### noduleMassWithH2assays.RDATA
##### h2withnods
	plant.strain
	plantID
	strain
	side.x
	nodsWeighed
	nodmassG
	phbnodFWg
	totalnodmassG
	nodcount
	file
	efficiencySlope
	r2
	EAC
	H2peak
	H2overCO2peak
	H2inair
	CO2inair
	CO2at0ppmH2
	date
	side.y
	excludeslope
	specificN2aseActivity
	identifier

##### h2withnods2
	strain
	plant.strain
	specificN2aseActivity
	H2peak
	identifier
	plantID
	side.x
	nodsWeighed
	nodmassG
	phbnodFWg
	totalnodmassG
	nodcount
	file
	efficiencySlope
	r2
	EAC
	H2overCO2peak
	H2inair
	CO2inair
	CO2at0ppmH2
	date
	side.y
	excludeslope
	rotation
	originalNodulePHBEst
	focal

##### h2withnodsw
	plantID
	strain.foc
	plant.strain.foc
	specificN2aseActivity.foc
	H2peak.foc
	identifier.foc
	side.x.foc
	nodsWeighed.foc
	nodmassG.foc
	phbnodFWg.foc
	totalnodmassG.foc
	nodcount.foc
	file.foc
	efficiencySlope.foc
	r2.foc
	EAC.foc
	H2overCO2peak.foc
	H2inair.foc
	CO2inair.foc
	CO2at0ppmH2.foc
	date.foc
	side.y.foc
	excludeslope.foc
	rotation.foc
	originalNodulePHBEst.foc
	strain.ref
	plant.strain.ref
	specificN2aseActivity.ref
	H2peak.ref
	identifier.ref
	side.x.ref
	nodsWeighed.ref
	nodmassG.ref
	phbnodFWg.ref
	totalnodmassG.ref
	nodcount.ref
	file.ref
	efficiencySlope.ref
	r2.ref
	EAC.ref
	H2overCO2peak.ref
	H2inair.ref
	CO2inair.ref
	CO2at0ppmH2.ref
	date.ref
	side.y.ref
	excludeslope.ref
	rotation.ref
	originalNodulePHBEst.ref
	specificN2aseActivity.normalized

##### nodMass
	strain
	plant.strain
	plantID
	side
	nodsWeighed
	nodmassG
	phbnodFWg
	totalnodmassG
	nodcount
	originalNodulePHBEst
	rotation
	masspernod
	focal

##### nodMassw
	plantID
	strain.foc
	plant.strain.foc
	side.foc
	nodsWeighed.foc
	nodmassG.foc
	phbnodFWg.foc
	totalnodmassG.foc
	nodcount.foc
	originalNodulePHBEst.foc
	rotation.foc
	masspernod.foc
	strain.ref
	plant.strain.ref
	side.ref
	nodsWeighed.ref
	nodmassG.ref
	phbnodFWg.ref
	totalnodmassG.ref
	nodcount.ref
	originalNodulePHBEst.ref
	rotation.ref
	masspernod.ref


#### phbplusplantdata_final.RDATA
##### mixedphbdata
	file
	FL1mean
	FL1cv
	FSCmean
	FSCcv
	SSCmean
	SSCcv
	countGated
	FL1mode
	FSCmode
	FL1geomean
	FSCgeomean
	pop
	runDate
	plate
	column
	row
	dilution
	note
	discard
	nodID

##### pd
	plot
	plantID
	tray
	experiment
	rotation
	yearswithoraftersoy
	sowDate
	harvestDate
	growthStage
	chl
	shootMassG
	nodMassG
	nodsWeighed
	nodCount
	note
	gpernodule
	nodMassG_plusPHBnods
	yearsofsoy
	tray1
	noduleMassCorrection
	correctedNodMassRatio
	oldornew
	rotation2

##### phbdata
	plantID
	plot
	file
	FL1mean
	FL1cv
	FSCmean
	FSCcv
	SSCmean
	SSCcv
	countGated
	countUngated
	FL1mode
	FSCmode
	FL1geomean
	FSCgeomean
	runDate
	plate
	column
	row
	dilution
	note
	discard
	phbslope
	phbintercept
	phbEstfrommean
	phbEst
	modephbEst
	tray
	experiment
	rotation
	yearswithoraftersoy
	sowDate
	harvestDate
	growthStage
	nodCount
	shootMassG
	nodID
	cellsPerNod
	yearsofsoy
	phbEst_zerobound
	plotrep
	plantrep
	oldornew
	rotation2
	percentileby10
	nodMassG
	gpernodule

##### standards
	sampleID
	file
	FL1mean
	FL1cv
	FSCmean
	FSCcv
	SSCmean
	SSCcv
	countGated
	countUngated
	FL1mode
	FSCmode
	FL1geomean
	FSCgeomean
	X
	runDate
	plate
	column
	row
	dilution
	note
	discard
	tubeNo
	phbPerCellPg
	cellsPerMl


##### uninoculatedplants
	plot
	plantID
	tray
	experiment
	rotation
	yearswithoraftersoy
	sowDate
	harvestDate
	growthStage
	chl
	shootMassG
	nodMassG
	nodsWeighed
	nodCount
	note
	gpernodule
	nodMassG_plusPHBnods

#### senescentNoduleDataEditedA3B1.RData
##### nr1
	nodID
	file
	sampleType
	stain
	dilutionFromTube
	note
	plantID
	strain
	FL3mean
	FSCmean
	FL3cv
	FSCcv
	Count
	X
	totalCellCount
	PIcellCount
	propDead
	phbEst

##### phist_expanded
	file
	FL3
	cellCount
	sampleType
	stain
	nodID
	dilutionFromTube
	note
	plantID
	strain
	phbEst

##### phist_expanded2
	file
	FL3
	cellCount
	sampleType
	stain
	nodID
	dilutionFromTube
	note
	plantID
	strain
	phbEst

##### standards
	file
	FL3mean
	FSCmean
	FL3cv
	FSCcv
	Count
	X
	sampleID
	sampleType
	stain
	dilutionFromTube
	note
	phbPerCellPg

##### standards_expanded
	file
	FL3
	cellCount
	FL3mean
	FSCmean
	FL3cv
	FSCcv
	Count
	X
	sampleID
	sampleType
	stain
	dilutionFromTube
	note
	phbPerCellPg
	phbEst
