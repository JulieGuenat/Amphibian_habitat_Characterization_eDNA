# Amphibian_habitat_Characterization_eDNA

## GENERAL INFORMATION


### 1. Title of Dataset 

Data from: Investigating fine scale breeding habitat use by amphibians in a continuous wetland using environmental DNA


### 2. Author Information

Julie Morgane Guenat1*, Antoine Gander2, Luca Fumagalli1,3§, Guillaume Lavanchy1,2§

1Laboratory for Conservation Biology, Department of Ecology and Evolution, University of Lausanne, 1015 Lausanne, Switzerland

2Association de la Grande Cariçaie, Chemin de la Cariçaie 3, 1400 Cheseaux-Noréaz, Switzerland

3Swiss Human Institute of Forensic Taphonomy, University Centre of Legal Medicine Lausanne-Geneva, Lausanne University Hospital and University of Lausanne, Ch. de la Vulliette 4, 1000 Lausanne 25, Switzerland

§joint last authors

*corresponding author: Julie.Guenat@ik.me


### 3. Date of data collection 

May 2018


### 4. Geographic location of data collection 

La Grande Cariçaie, Lake Neuchâtel, Switzerland


## DATA & FILE OVERVIEW


### 1. Description of dataset

We used eDNA metabarcoding to study the fine scale breeding habitat use of amphibians in a continuous wetland expanse in Switzerland. We characterize the fine scale breeding habitat of the seven amphibian species present there (Lissotriton vulgaris, L. helveticus, Bufo bufo, Rana temporaria, Hyla arborea and Pelophylax ridibundus, P. bergeri) by testing the influence of six abiotic environmental variables (average mud depth, percentage of emerged and submerged vegetation, percentage of emerged land, average water temperature and the distance to the wintering habitat) on community structure and individual species distribution. Finally, we investigated whether biotic interactions could influence amphibian distribution by examining co-occurrences of the different species in the studied area.


Sample collection:

We collected water samples from May 21st to 28th 2018, which corresponds to the breeding season for local amphibians. We collected two liters of water at each of the 50 sampling points (i.e., 25 in Yverdon and 25 in Gletterens) using the VigiDNA kit (Spygen, Le Bourget du Lac, France) following manufacturer's protocol. Filtration capsules were stored at room temperature for two months, until we extracted DNA.


### 2. File List: 
	
File name: OBItools_pipeline.sh

Description: script obitools (commands and parameter) used to process raw barcode sequences. 


File name: 12S_P.bergeri_seq_682bp.txt

Description: consensus partial sequence of the 12S mitochondrial gene of Pelophylax bergeri. Amplicon amplified with primers L2519 and H3296 (Wang et al. 2017) was then sequenced using Sanger sequencing. 


File name: Statistical_pipeline.Rmd

Description: Rmarkdown describing the statistical analyses performed. 


File name: final_env_Dec20.csv

Description: Environmental data collected at each sampling points. The data were collected from April 21st to 23rd 2018. 


File name: 8_Batr_presabs.txt

Description: Infered presence/absence data of the six amphibian species at each sampling point.  
