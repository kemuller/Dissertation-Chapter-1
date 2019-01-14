## Analysis scripts used for manuscript:

### Title: Resource acquisition and allocation traits in symbiotic rhizobia with implications for life-history outside of legume hosts

#### Authors:
Katherine E. Muller (collected data)
R. Ford Denison

#### Date of data collection: 2015 - 2018

#### Geographic location of data collection (where was data collected?): St. Paul, MN

#### Information about funding sources that supported the collection of the data:
Both authors were funded by the Minnesota Long-Term Agricultural Research Network. KEM was also funded by the Carol Pletcher Fellowship and Carolyn Crosby Grant (University of Minnesota Graduate School). 

--------------------------
Script:
## PHBvariationInNaturalPopulations_edited.R

Requirements:
packages: plyr, ggplot2, lme4, sjstats, lattice
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
Measurements of Ensifer meliloti standards used to estimate PHB/cell from fluorescence of a lipid stain measured in the flow cytometer. This goes with flowDataA3B5nod_bdpy.csv
	file: name of raw data file obtained from the flow cytometer
	FLmean: geometric mean of PHB fluorescence from the flow cytometry data (stained with Bodipy, calculated in FlowJo).
	gatedCells: number of cells included in FLmean (gated by hand in FlowJo to exclude debris)
	date: date samples were measured on the flow cytometer
	stain: stain used for standards. levels are "nile red" and "bodipy" (names of different stains used to measure PHB).
	standardID: ID for the E. meliloti standard
	PHBperCellpg: PHB per cell measured previously with gas chromatography by Will Ratcliff.
	standardID2: I originally made a mistake when entering names for the standards. this column is corrected by hand. 
	PHBperCellpg2: corresponds to corrected standard ID.

#### flowDataA3B5nod_bdpy.csv
	file: name of raw data file obtained from the flow cytometer. Each file corresponds to a single root nodule.
	FL1mean: geometric mean of PHB fluorescence from the flow cytometry data (stained with Bodipy, calculated in FlowJo).
	gatedCells: number of cells included in FL1mean (gated by hand in FlowJo to exclude debris)
	date: date samples were measured on the flow cytometer
	plantID: ID for the plant nodules were harvested from (soybean cultivar MN0095 in growth pouches).
	nodAgeDays: days between when the nodule was formed (recorded on growth pouch) and when it was harvested)
	nodLet: letter assigned nodules on each plant (combined with plantID to make unique nodule ID). 
	dilution: how much stained sample was diluted during flow cytometer (100 = 100-fold dilution, 10 = 10-fold dilution, 1 = no dilution).
	note: notes made during flow cytometry
	delete: I used two dilutions for some samples. I marked a 1 in this column to indicate the sample I decided not to use because it was either too concentrated or too dilute.

#### averagemonthlysoiltemp_waseca.csv
Soil temperature data from the Southern Research and Outreach Center in Waseca, MN. Entered from historic weather reports: https://sroc.cfans.umn.edu/weather-sroc/historic-reports. Soil temperature is measured daily (max and min). In the model, I used their reported average soil temperature for each month and year measured at a depth of 8 inches (slightly more stable than soil temperature at 4 inches). They reported it in Farenheight and I converted it to Celsius.
	year
	month
	average4inchtempF
	average8inchtempF
	average4inchtempC
	average8inchtempC

#### contaminationNST1.csv
Data on contaminants observed in starvation cultures.
	sample: The number corresponds to consecutive samples over time (1st, 2nd, etc.). The order of the sample # corresponds to the date.
	date: date on which aliquots of starvation culture were sampled (plated for viability counts and fixed in ethanol for PHB).
	density: LD = low-density, HD = high-density (starvation cultures ~10^5 cells/mL or ~10^5 cells/mL based on OD600)
	strain: code used to identify isolate in starvation culture (110 = inoculum strain Bradyrhizobium diazoefficiens USDA110. A3, A4, etc. are isolates from soybean nodules from Waseca, MN). 
	type: morphotypes of contaminants. Codes I, II, II, IV indicate different morphotypes. I gave up on morphotyping after the first few samples.
	count1: count of colonies within a plate or drop 
	r1: indicates dilution of the sample used for count1 (row of the 96-well plate used for 5-fold dilution. Usually 1 = 5^0, 2 = 5^1, etc.). 
	count2: count of colonies within a second plate or drop (if available)
	r2: dilution used for count2 (see above)
	count3: count of colonies within a third plate or drop (if available). 
	r3: dilution used for count3 (see above)
	v1: volume of sample used for count1 (in uL, used to estimate colonies/mL)
	v2: same as v1, but for count2
	v3: same as v1, but for count3
	rep: replicate aliquot (I usually sampled two aliquots at a time).



#### isolateCodesForFreezerwithfilenames_WRE.csv
This provides information about single-colony isolates from Waseca used in split-root plants. 
	nodID: ID of the nodule the isolate came from (matches nodID in phbplusplantdata_final.RDATA
	file: file from flow cytometer
	phbEst: PHB/cell (pg) estimated from S. meliloti standards
	code: alphanumeric code used to identify isolate
	repFromGrowthCurves: I measured growth curves in tubes while I was growing single-colony isolates (OD6000). That required replicate tubes. I only archived one of the replicates (the number doesn't really matter).
	plantID: ID of the plant the isolate came from (from phbplusplantdata_final.RDATA)
	treatment: describes history of soil where the isolate came from (30-year soybean or 30-year corn plots in Waseca).
	rotation: shortened codes for treatment column
	phbEstfromArithmeticMean: I originally estimated PHB from arithmetic mean PHB fluorescence (instead of geometric mean) by mistake. I put the corrected estimate in the phbEst column and kept the column with the arithmetic mean so I would know that the mistake had been fixed.

#### chamecristaCombinedNodData_edited.RDAtA:
##### nod.data.combined
	plantID: ID of Chamaecrista fasciculata plant from the field(GA01 = greycloud dunes, transect A, 1 meter mark, GB022 = greycloud dunes, transect B, meter 22). I originally had two sites, but none of my plants survived on the other site. 
	nodID: unique identifier for each nodule (plant ID plus an integer). 
	file: raw data file from flow cytometer (edited in FlowJo--measurements below are from flow jo, for gated cells only).
	FL1mean: geometric mean PHB fluorescence (Bodipy measured on FL1 channel)
	FL1cv: coefficient of variation for FL1
	FSCmean: geometric mean forward scatter (proxy of cell size)
	FSCcv: coefficient of variation for forward scatter
	SSCmean: geometric mean side-scatter (proxy of cell shape)
	SSCcv: coefficient of variation for side-scatter
	countGated: count of cells gated by hand to exclude debris (in flowjo). 
	FL1mean_high: some nodule samples had distinguishable high- and low-PHB groups (perhaps mixed infection?). These were gated by hand in flowjo. The names with _high and _low are the same as above, except for gated groups. 
	FSCmean_high
	SSCmean_high
	FL1mean_low
	FSCmean_low
	SSCmean_low
	countGated_high
	countGated_low
	runDate: date samples were measured on the flow cytometer.
	dilutionFromNodule: 100 means fixed samples were diluted 100-fold from original nodule extract.
	dilutionFromTube: 1 means fixed, stained samples were not diluted when measured on the flow cytometer
	flowRate: flow rate setting on the flow cytometer (low or med). There were some technical problems with the low flow rate setting. 
	note: notes made during flow cytometry
	phbEst: estimated PHB/cell from FL1mean
	wt.g: fresh weight of nodules (in g)
	cfuperml: estimated CFU/mL from plate counts on AG agar
	propSmooth: originally I thought there were different colony morphotypes(smooth vs. less smooth), but it turned out to be changes in colony appearance based on age.
	growsWithBG: indicates whether or not the nodule extract grew in AG agar supplemented with 1ug/mL Brilliant Green (dilution plates were done with and without BG). 
	growthStage: describes maturity of plants when they were harvested (flw = flowers present, earlypod = some pods present, not fully formed; bud = flower buds present, but not opened; veg = no buds or flowers).
	shootDWg: mass of dried shoots (in grams).

#### editedPHBData_WREsplitroot1.RDATA
PHB measurements from splitroot plants (Fig. 1 in paper). All flow cytometry measurements (FL1, FSC, SSC) are for cells hand-gated to exclude debris (using flowjo). FL1 = PHB flourescence (Bodipy measured on the FL1 channel). FSC = forward scatter, SSC = side scatter. mean = arithmetic mean for thousands of cells. Mode = mode for thousands of cells. geomean = geometric mean for thousands of cells (this was what I actually used for estimating PHB/cell). 
##### phbdata
	file: raw data file from flow cytometer (each file has measurements for thousands of cells from a single nodule)
	FL1geomean	FL1mode
	FL1cv
	FSCmean
	FSCmode
	FSCcv
	SSCmean
	SSCmode
	SSCcv
	countGated: count of cells gated to exclude debris
	countUngated: count of all cells in sample (including debris)
	FSCgeomean
	FL1mean
	X
	runDate
	plate: measurement block (96-well plate in which samples were fixed, stained, and measured). The plate and position on the plate (row and column) were used to keep track of samples during processing (also acted as a blind).
	row: row on 96-well plate
	column: column on 96-well plate. 
	strain: there should be two strains per plant. The strain ID correspond to the isolate code in isolateCodesForFreezerwithfilenames_WRE.csv. The ID for the common reference strain is A3. 
	plantID: indicates split-root plant nodules were harvested from.
	skip: I tried some samples at more than one dilutions or flow rates. SKIP means I decided not to use that sample
	dilution: dilution factor from stained sample (1 = no dilution, 10 = 10-fold dilution, etc.)
	flowRate_ulpermin: flow rate setting on the flow cytometer (12 = low, 35 = med). 
	note
	phbEst: PHB/cell estimated from FL1geomean and E. meliloti standards (in dataframe standards)


##### standards
Flow cytometry measurements for E. meliloti standards measured alongside samples in phbdata. Columns are the same as phbdata except where indicated. 
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
	phbPerCellPg: PHB/cell (pg) measured with gas chromatography by Will Ratcliff. 


#### noduleMassWithH2assays.RDATA
These are additional measurements from the splitroot plants (Fig 1C). Some of the measurements are from a gas exchange assay that isn't discussed in this paper. At the time it seemed easier to keep them in than to try and pick out only the columns I needed to use.  
##### h2withnods
Ignore this dataframe (it's a slightly rougher version of h2withnods2)

##### h2withnods2
long-format version of h2withnodsw (each plant has 2 rows, one for each side of the pouch)
	strain: ID of isolate (A3 = reference strain; A1, B1, etc. = focal isolates)
	plant.strain: unique idenfier for the side of the plant (long format)
	specificN2aseActivity: from gas exchange (ignore)
	H2peak: from gas exchange (ignore)
	identifier: unique ID for nodule
	plantID: ID of split-root plant nodules were harvested from
	side.x: from gas exchange (ignore)
	nodsWeighed: number of nodules weighed after plants were dried
	nodmassG: mass of dried nodules (g, excludes nodules harvested for PHB)
	phbnodFWg: fresh mass of nodules harvested for PHB (g). 
	totalnodmassG: mass of dried nodules (in g) plus estimated dried mass of nodules harvested from PHB (fresh weight * 0.285--based on lab measurements of fresh and dried soybean nodules)
	nodcount: nodules per plant
	file: from gas exchange (ignore)
	efficiencySlope: from gas exchange (ignore)
	r2: from gas exchange (ignore)
	EAC: from gas exchange (ignore)
	H2overCO2peak: from gas exchange (ignore)
	H2inair: from gas exchange (ignore)
	CO2inair: from gas exchange (ignore)
	CO2at0ppmH2: from gas exchange (ignore)
	date: from gas exchange (ignore)
	side.y: from gas exchange (ignore)
	excludeslope: from gas exchange (ignore)
	rotation
	originalNodulePHBEst: PHB/cell in source nodule for that isolate (see isolateCodesForFreezerwithfilenames_WRE.csv)
	focal: ref = reference strain (A3), foc = focal strain (A1, B1, C1, etc.)

##### h2withnodsw
same as h2withnods2 except in wide format (same columns) Each plant has one row with separate columns for each side of the pouch (foc = focal isolate, ref = common reference strain). 
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
	specificN2aseActivity.normalized: ignore (gas exchange)

##### nodMass
nodule mass data from splitroot plants (same as columns described in h2withnods2, except where indicated). Long-format version (each plant has two rows)
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
	masspernod: nodules per plant divided by total nodule mass per plant.
	focal

##### nodMassw
same as nodMass except in wide format (each plant has one row)
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
PHB measurements from trap-plants inoculated with field soil (Figure 1B). One row per nodule.
##### phbdata
	plantID: ID of trap plant (from greenhouse)
	plot: ID of field plot where soil was collected in Waseca, MN
	file: name of data file from the flow cytometer (a file = a nodule)  
Flow cytometry summary stats from Flowjo (applied to subsets gated to exclude debris). FL1 = PHB fluorescence (Bodipy measured on the FL1 channel), FSC = forward scatter, SSC = side scatter. mean = arithmetic mean, mode = mode (value with highest density), geomean = geometric mean. I was unsure at the time which summary stat to use, so I decided to collect all of them to compare. cv = coefficient of variation. 
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
	FL1geomean: geometric mean PHB fluorescence (used to estimate PHB/cell)
	FSCgeomean
	runDate: date samples were measured on the flow cytometer
	plate: measurement block (96-well plate used to keep track of samples). 
	column: horizontal position within measurement block
	row: vertical position within measurement block
	dilution: how much the stained sample was diluted before flow cytometry (100 = 100-fold, 10 = 10-fold, etc.). 
	note: notes added during flow cytometry
	discard: NA 
	phbslope: slope of standard curve used to estimate PHB/cell from E. meliloti standards (separate set of standards run for each plate)
	phbintercept: intercept of standard curve used to estimate PHB/cell from E. meliloti standards
	phbEstfrommean: PHB estimate from arithmetic mean (ignore)
	phbEst: PHB/cell (pg) estimated from FL1geomean with phbslope and phbintercept.
	modephbEst: PHB/cell estimated with the mode (ignore)
	tray: indicates blocks of plants in the greenhouse.
	experiment: codes indicating which field experiment the soil came from (two different field experiments in Waseca, doesn't matter for the sake of this paper). 
	rotation: codes for cropping history of field plots (not a part of this paper)
	yearswithoraftersoy: rotation coded differently (ignore)
	sowDate: date trap plants were sown (used the date dry seeds became wet)
	harvestDate: date trap plants were harvested. 
	growthStage: 2 = pod-filling, 1 = flowering, 0 = not-yet flowering.
	nodCount: nodules per plant
	shootMassG: mass of dried shoots (in g)
	nodID: individual identifier for the nodule (matches with isolateCodesForFreezer)
	cellsPerNod: cell count approximated from flow cytometry data (not reliable--I learned later that there were issues with the flow rate).
	yearsofsoy: yet another re-coding of rotation
	phbEst_zerobound: same as PHBest, except negative values are set to zero (used for Fig. 1B).
	plotrep: field plots coded as 1-4 (ignore)
	plantrep: replicate trap plants (inoculated with soil from the same plot) coded as 1-4 (ignore)
	oldornew: another way of categorizing rotation treatment (ignore)
	rotation2: yet another re-coding of rotation (ignore)
	percentileby10: I don't remember what this is. But it's not something I use or talk about in the paper (ignore).
	nodMassG: dry mass of nodules/plant (excluding nodules harvested for PHB).
	gpernodule: nodMassG/nodCount. 

##### standards
E. meliloti standards used to estimate PHB/cell. The columns are mostly the same as phbdata. They match up with the "plate" column (a set of standards was stained and run alongside each plate). 
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
	phbPerCellPg: PHB/cell (pg) measured with gas chromatography by Will Ratcliff. 
	cellsPerMl

##### mixedphbdata
There were 107 nodules that looked like they had distinct subpopulations differing in PHB/cell. I gated the subpopulations in Flowjo and calulated separate summary statistics for each. It was a fairly small portion of the 1276 nodules presented in Fig. 1B, so I didn't go into it in the paper. The columns are the same except where indicated. It's in a long format, so each subpopulation within a nodule gets its own row (indicated in 'pop').
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
	pop: high, med, low = distinct subpopulations with high, medium, or low PHB (gated with flowjo).  
	runDate
	plate
	column
	row
	dilution
	note
	discard
	nodID

##### pd
Data on trap plants inoculated with field soil (one row per plant)
	plot: ID of field plot where soil was collected in Waseca, MN
	plantID: unique identifier for the plant (includes plot name)
	tray: randomized block in the greenhouse (12 plants were grouped onto a tray). 
	experiment: LTARN and SROC indicate different field experiments in Waseca MN (doesn't matter for this paper).
	rotation: code for cropping history (doesn't matter for this paper)
	yearswithoraftersoy: another way of coding cropping history (ignore)
	sowDate: date plants were sown in the green house (I use date seeds went from dry to wet)
	harvestDate: date plants were harvested (nodules collected)
	growthStage: 2 = pod-filling, 1 = flowering, 0 = not-yet flowering.
	chl: leaf measurements with AtLeaf Chlorophyll meter (ignore)
	shootMassG: dry mass of shoots (g)
	nodMassG: dry mass of nodules in grams (not including nodules harvested for PHB)
	nodsWeighed: nodules weighed for dry mass
	nodCount: nodules per plant
	note: notes taken while weighing plant biomass
	gpernodule: nodMassG divded by nodsWeighed
	nodMassG_plusPHBnods: nodule mass with estimated DW of nodules harvested for PHB added back in (gpernodule * 8). I ended up not using this. 
	yearsofsoy: another way of coding cropping history (ignore)
	tray1: originally the idea was to have no replicates from the same plot on the same tray and to use the plotID and tray number to make the plantID. But of course, mistakes were made. This column was used instead.
	noduleMassCorrection: this was the estimated DW of nodules harvested for PHB.
	correctedNodMassRatio: I don't remember what this was, but it's not used for anything
	oldornew: another way of coding cropping history 
	rotation2: another way of coding cropping history. 



##### uninoculatedplants
Same as pd, except for uninoculated controls in the greenhouse. None of them had nodules. 
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
PHB measurements from senescent soybean nodules used to make histogram in Figure 2A. 
##### nr1
	nodID: unique identifier for the nodule
	file: name of data file from flow cytometer (1 per nodule sample)
	sampleType: indicates they are senescent nodules (original data file also included E. meliloti standards, which are separated into a different dataframe)
	stain: indicates they were stained with Nile Red
	dilutionFromTube: 1 means fixed, stained samples were measured at full strength on the flow cytometer (no additional dilution)
	note: notes made during flow cytometry
	plantID: ID of plant nodules were harvested from.
	strain: identity of isolate (A3 is the reference strain from split root plants, B1 is an isolate from Waseca used as one of the focal isolates in split-root plants). 
	FL3mean: geometric mean PHB fluorescence (from Nile Red measured on the FL3 channel--stats for gated cells only).
	FSCmean: geometric mean forward scatter
	FL3cv: coefficient of variation for PHB fluorescence
	FSCcv: coefficient of variation for forward scatter
	Count: count of gated cells (debris removed by hand in flowjo).
	X: stupid blank column that Flowjo inserts sometimes
	totalCellCount: count of gated (non-debris) cells in live samples stained with propidium iodide (dead cells absorb PI, live ones do not)
	PIcellCount: subset of totalCellCount that appeared to absorb propidium iodide (measured on the FL3 channel). Gated in flowjo. This didn't work very well.
	propDead: PIcellCount divided by totalCellCount
	phbEst: estimated PHB per cell from FL3mean (estimated with E. meliloti standards). 

##### phist_expanded
Fluoresecence data from individual cells (exported from histograms in Flowjo). Except where indicated, column names match those described in nr1. 
	file
	FL3: Nile red fluorescence measured on the FL3 channel
	cellCount: count of cells with the specified FL3 value.
	sampleType
	stain
	nodID
	dilutionFromTube
	note
	plantID
	strain
	phbEst: PHB estimated from FL3 using E. meliloti standards.

##### phist_expanded2
Same as phist_expanded except it removes one file with bad data) (NRnodules.016)
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
Ensifer meliloti standards used to estimate mean PHB/cell in nodule samples (columns are the same as in nr1 except where indicated). 
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
	phbPerCellPg: PHB/cell (pg) measured with gas chromatography by Will Ratcliff.

##### standards_expanded
Same as phist_expanded except for the E. meliloti standards
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
